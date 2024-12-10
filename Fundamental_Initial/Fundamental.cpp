// Imagine++ project
// Project:  Fundamental


//kharef okba

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}
// Method to select a random subset of 8 matches
vector<int> selectRandomSubset(vector<Match>& matches) {
    vector<int> inds;
    // Seed the random number generator with the current time
    //srand(static_cast<unsigned int>(time(0)));
    // Select 8 random indices
    for (int i = 0; i < 8; ++i) {
        int temp = rand() % matches.size();
        inds.push_back(temp);
    }
    return inds;
}


FMatrix<float,3,3> tempF(vector<Match>& matches, vector<int>& inds2 ,bool ref){
    Matrix<float> A;
    float x1, y1, x2, y2;
    float No= 0.001;
    vector<int> inds;
    if(ref){
        inds.assign(inds2.begin(), inds2.end());
        A=Matrix<float>(inds2.size(),9);
        for (int i=0; i<inds.size(); i++){
            //cout<<i<< " :index for refining "<<endl;
            // Normalization
            x1 = No * matches[inds[i]].x1;
            y1 = No * matches[inds[i]].y1;
            x2 = No * matches[inds[i]].x2;
            y2 = No * matches[inds[i]].y2;
            // Building the linear system
            A(i,0) = x1 * x2;
            A(i,1) = x1 * y2;
            A(i,2) = x1;
            A(i,3) = y1 * x2;
            A(i,4) = y1 * y2;
            A(i,5) = y1;
            A(i,6) = x2;
            A(i,7) = y2;
            A(i,8) = 1;
        }
         A(8,8) = 0;
    }
    else{
         A=Matrix<float>(9,9);
         inds= selectRandomSubset(matches);
         for (int i=0; i<inds.size(); i++){
            //cout<<i<< " index random subset "<<endl;
            // Normalization
            x1 = No * matches[inds[i]].x1;
            y1 = No * matches[inds[i]].y1;
            x2 = No * matches[inds[i]].x2;
            y2 = No * matches[inds[i]].y2;
            // Building the linear system
            A(i,0) = x1 * x2;
            A(i,1) = x1 * y2;
            A(i,2) = x1;
            A(i,3) = y1 * x2;
            A(i,4) = y1 * y2;
            A(i,5) = y1;
            A(i,6) = x2;
            A(i,7) = y2;
            A(i,8) = 1;
         }
         A(8,8) = 0;

    }
    // Solve linear system using svd
    Vector<float> S;
    Matrix<float> U, V;
    svd(A,U,S,V);

    FMatrix<float,3,3> F;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            F(i,j) = V(8,3*i+j);
        }
    }

    FVector<float,3> S2;
    FMatrix<float,3,3> U2, V2;
    svd(F,U2,S2,V2);

    S2[2] = 0;

    FMatrix<float,3,3> N(0.f);
    N(0,0) = 0.001;
    N(1,1) = 0.001;
    N(2,2) = 1;

    F = transpose(N)* U2 * Diagonal(S2) * V2 * N;

    return F;
}

float CalcDistance(float x1,float y1,float x2,float y2,FMatrix<float,3,3> F){
    FVector<float,3> x   = FVector<float,3>(x1,y1,1);
    FVector<float,3> x0   = FVector<float,3>(x2,y2,1);

    FVector<float,3> point = transpose(F) * x;
    float dist = abs(x0*point);
    dist= dist / sqrt(point[0]*point[0] + point[1]*point[1]);
    return dist;

}

// Given F and the matches, returns the indices of inliers for F
vector<int> findInliers(FMatrix<float,3,3>& F, vector<Match>& matches, float distMax){
    vector<int> inliers;
    float x1, x2, y1, y2;
    for (int i = 0; i < matches.size(); i++){
        x1 = matches[i].x1;
        y1 = matches[i].y1;
        x2 = matches[i].x2;
        y2 = matches[i].y2;


        float dist=CalcDistance(x1,y1,x2,y2,F);
        if (dist <= distMax){
            inliers.push_back(i);
        }
    }
    return inliers;
}


// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    vector<int> bestInliers;
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF POINTS
    vector<int> Inliers ;
    int nInliersMax = 0;

    int n = matches.size();
    for(int iter=0; iter<Niter; iter++) {

        FMatrix<float,3,3> F = tempF(matches, Inliers,false);
        Inliers = findInliers(F, matches, distMax);

        if (Inliers.size() > nInliersMax) {
            nInliersMax = Inliers.size();
            bestF = F;
            bestInliers = Inliers;

            // Update estimate of Niter
            float InliersOP = Inliers.size() / static_cast<float>(n);
            // Safeguard against invalid values for InliersOP
            // Ensure InliersOP is within the valid range (0 < InliersOP < 1) for log and pow calculations
            if (InliersOP <= 0.0f) {
                InliersOP = std::numeric_limits<float>::min();  // Set to a very small positive value
            } else if (InliersOP >= 1.0f) {
               InliersOP = 1.0f - std::numeric_limits<float>::epsilon();  // Slightly less than 1 to avoid log(0)
            }
            // Safeguard against invalid values for log and pow calculations
            if (InliersOP > 0 && InliersOP < 1) {
                // Protect against division by zero and ensure valid log inputs
                Niter = min(Niter, static_cast<int>(log(BETA) / log(1 - pow(InliersOP, 8))));
            }
        }
        cout<<iter<< " :Number of iterations "<<endl;
    }


    std::cout << "Size of bestInliers: " << bestInliers.size() << std::endl;
    std::cout << "Size of matches: " << matches.size() << std::endl;

    //refining F
    bestF =tempF(matches,bestInliers,true);
    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    return bestF;
}

// Draw a cross at position (x,y)
void drawCross(int x, int y, Color color, int size = 5, int penWidth=1){
    drawLine(x-size, y-size, x+size, y+size, color, penWidth);
    drawLine(x-size, y+size, x+size, y-size, color, penWidth);
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {

    while(true) {
        int x,y;
        FVector<float, 3>  p1, l;
        int w1 = I1.width();
        int w2 = I2.width();
        if(getMouse(x,y) == 3)
            break;
        Color color=YELLOW ;
        drawCross(x,y,color);

        if (x < w1){ // Click in left image
            color=RED;
            p1 = FVector<float, 3>(x, y, 1.0);  // Point clicked in first image
            l = transpose(F) * p1;
            float y = -l[2]/l[1];
            float y0 = -(l[0]*(float)(w2) + l[2])/l[1];
            drawLine(w1,(int)y,w1+w2,(int)y0,color);
        }
        else { // Click in right image
            color=BLUE;
            p1 = FVector<float, 3>((x - w1), y, 1.0); //point clicked in the second image
            l = F * p1;
            float y0 = -l[2]/l[1];
            float y1 = -(l[0]*(float)w1 + l[2])/l[1];
            drawLine(0,(int)y0,w1,(int)y1,color);
        }
    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    cout << " matches: " << matches.size() << endl;
    click();

    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);
    }
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
