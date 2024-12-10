// Imagine++ project
// Project:  Panorama
// okba kharef

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <sstream>
using namespace Imagine;
using namespace std;

// Function to draw a cross at the specified (x, y) coordinates.
// The cross is drawn using two diagonal lines centered at (x, y).
// - 'size' defines half the length of the lines.
// - 'penWidth' is the thickness of the lines.
void drawCross(int x, int y, Color color, int size = 5, int penWidth = 2) {
    // Draw the first diagonal line from top-left to bottom-right
    drawLine(x - size, y - size, x + size, y + size, color, penWidth);

    // Draw the second diagonal line from bottom-left to top-right
    drawLine(x - size, y + size, x + size, y - size, color, penWidth);
}

// Function to record clicks on two images and match pairs of points.
// It runs until the right mouse button is clicked.
void getClicks(Window w1, Window w2,
               vector<IntPoint2>& pts1, vector<IntPoint2>& pts2) {
    // ------------- TODO/A completer ----------
    // Create a guide window to display instructions and the point-pair count.
    Window Guide = openWindow(950, 200);  // Open a window with a size of 950x200 pixels
    setActiveWindow(Guide);  // Set the Guide window as the active window to draw on

    // Display initial instructions for the user in the Guide window
    drawString(10, 30, "Select at least 4 pairs of matching points in the two images.", BLACK);
    drawString(10, 65, "Pairs count :", BLACK);
    drawString(280, 65, "0", BLACK);  // Initialize the pair count display to 0

    // Variables to track if a click has occurred in each window
    bool clickedW1, clickedW2 = false;
    int i = 0;
    IntPoint2 click_coordinates;
    int sw;
    Window W;

    // Main loop to capture clicks until the right mouse button is clicked
    while (true) {
        // Capture a mouse click event. The clicked window is stored in W,
        // and the click coordinates are stored in 'click_coordinates'.
        int button = anyGetMouse(click_coordinates, W, sw);

        // If the right mouse button (button == 3) is clicked, break the loop and stop collecting points.
        if (button == 3) {
            break;
        } else {
            // If the click was in the first window (w1)
            if (W == w1) {
                pts1.push_back(click_coordinates);  // Add the point to the 'pts1' vector
                setActiveWindow(w1);  // Set window 1 as the active window
                drawCross(click_coordinates.x(), click_coordinates.y(), RED);  // Draw a red cross at the clicked point
                clickedW1 = true;  // Mark that window 1 has received a click
            }
            // If the click was in the second window (w2)
            else if (W == w2) {
                pts2.push_back(click_coordinates);  // Add the point to the 'pts2' vector
                setActiveWindow(w2);  // Set window 2 as the active window
                drawCross(click_coordinates.x(), click_coordinates.y(), RED);  // Draw a red cross at the clicked point
                clickedW2 = true;  // Mark that window 2 has received a click
            }

            // If both windows (w1 and w2) have received clicks, we have a valid pair of points
            if (clickedW1 && clickedW2) {
                i++;  // Increment the pair count
                clickedW1 = false;  // Reset the click flag for window 1
                clickedW2 = false;  // Reset the click flag for window 2

                // Update the Guide window with the new pair count immediately
                setActiveWindow(Guide);  // Set the Guide window as the active window
                fillRect(250, 35, 315, 80, WHITE);  // Clear the area where the number is displayed
                drawString(280, 65, to_string(i), BLACK);  // Display the updated pair count
            }
        }
    }
}


// Return homography compatible with point matches
Matrix<float> getHomography(const vector<IntPoint2>& pts1,
                            const vector<IntPoint2>& pts2) {
    size_t n = min(pts1.size(), pts2.size());
    if(n<4) {
        cout << "Not enough correspondences: " << n << endl;
        return Matrix<float>::Identity(3);
    }
    Matrix<double> A(2*n,8);
    Vector<double> B(2*n);
    // ------------- TODO/A completer ----------

    for (size_t j=0 ; j<n ; j++){
        // Retrieve two corresponding points (x1,y1) and (x2,y2)
        double x1 = pts1[j].x();
        double y1 = pts1[j].y();
        double x2 = pts2[j].x();
        double y2=  pts2[j].y();

        //fill in matrix A and  vector B based on retireved points
        B[2*j]=x2;
        B[2*j+1]=y2;

        A(2*j, 0)= x1;
        A(2*j, 1)= y1;
        A(2*j, 2)= 1;
        A(2*j, 3)= 0;
        A(2*j, 4)= 0;
        A(2*j, 5)= 0;
        A(2*j, 6)= -x2*x1;
        A(2*j, 7)= -x2*y1;

        A(2*j+1, 0)= 0;
        A(2*j+1, 1)= 0;
        A(2*j+1, 2)= 0;
        A(2*j+1, 3)= x1;
        A(2*j+1, 4)= y1;
        A(2*j+1, 5)= 1;
        A(2*j+1, 6)= -y2*x1;
        A(2*j+1, 7)= -y2*y1;
    }
    B = linSolve(A, B);
    Matrix<float> H(3, 3);
    H(0,0)=B[0]; H(0,1)=B[1]; H(0,2)=B[2];
    H(1,0)=B[3]; H(1,1)=B[4]; H(1,2)=B[5];
    H(2,0)=B[6]; H(2,1)=B[7]; H(2,2)=1;

    // Sanity check
    for(size_t i=0; i<n; i++) {
        float v1[]={(float)pts1[i].x(), (float)pts1[i].y(), 1.0f};
        float v2[]={(float)pts2[i].x(), (float)pts2[i].y(), 1.0f};
        Vector<float> x1(v1,3);
        Vector<float> x2(v2,3);
        x1 = H*x1;
        cout << x1[1]*x2[2]-x1[2]*x2[1] << ' '
             << x1[2]*x2[0]-x1[0]*x2[2] << ' '
             << x1[0]*x2[1]-x1[1]*x2[0] << endl;
    }
    return H;
}

// Grow rectangle of corners (x0,y0) and (x1,y1) to include (x,y)
void growTo(float& x0, float& y0, float& x1, float& y1, float x, float y) {
    if(x<x0) x0=x;
    if(x>x1) x1=x;
    if(y<y0) y0=y;
    if(y>y1) y1=y;
}

// ------------- HELPER METHODS: Wrap and blend the two images ----------

// Helper method to transform coordinates using the inverse homography
Vector<float> transformCoordinates(const Vector<float>& p, const Matrix<float>& Hinv) {
    Vector<float> point_transformed = Hinv * p;
    point_transformed /= point_transformed[2];  // Normalize homogeneous coordinates
    return point_transformed;
}

// Helper method to interpolate the color from the image at given coordinates
Color interpolateImage(const Image<Color, 2>& img, const Vector<float>& coords) {
    return img.interpolate(coords[0], coords[1]);
}

// Helper method to blend two colors
Color blendColors(const Color& color1, const Color& color2) {
    Color blended;
    blended.r() = (color1.r() + color2.r()) * 0.5;
    blended.g() = (color1.g() + color2.g()) * 0.5;
    blended.b() = (color1.b() + color2.b()) * 0.5;
    return blended;
}

void blendImages(const Image<Color, 2>& I1, const Image<Color, 2>& I2,
                 Image<Color>& panorama, const Matrix<float>& Hinv,
                 float x0, float y0) {
    Vector<float> point(3);

    for (int y = 0; y < panorama.height(); y++) {
        for (int x = 0; x < panorama.width(); x++) {
            // Calculate the panorama coordinates
            point[0] = x0 + x;
            point[1] = y0 + y;
            point[2] = 1;

            // Check if pixel is within I2 bounds and interpolate
            if (point[0] >= 0 && point[0] < I2.width() && point[1] >= 0 && point[1] < I2.height()) {
                panorama(x, y) = interpolateImage(I2, point);
            }

            // Transform panorama coordinates to I1 using the inverse homography
            Vector<float> point_transformed = transformCoordinates(point, Hinv);

            // Check if transformed coordinates are within I1 bounds
            if (point_transformed[0] >= 0 && point_transformed[0] < I1.width() &&
                point_transformed[1] >= 0 && point_transformed[1] < I1.height()) {

                Color I1_color = interpolateImage(I1, point_transformed);

                // Blend the colors if the pixel is also within I2 (overlapping region)
                if (point[0] >= 0 && point[0] < I2.width() && point[1] >= 0 && point[1] < I2.height()) {
                    panorama(x, y) = blendColors(panorama(x, y), I1_color);
                } else {
                    // Use I1's color if not overlapping
                    panorama(x, y) = I1_color;
                }
            }
        }
    }
}

// Panorama construction
void panorama(const Image<Color,2>& I1, const Image<Color,2>& I2,
              Matrix<float> H) {
    Vector<float> v(3);
    float x0=0, y0=0, x1=I2.width(), y1=I2.height();

    v[0]=0; v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=0; v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    cout << "x0 x1 y0 y1=" << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1<<endl;

    Image<Color> I(int(x1-x0), int(y1-y0));
    setActiveWindow( openWindow(I.width(), I.height()) );
    I.fill(WHITE);
    // ------------- TODO/A completer ----------

    // Calculate the inverse of homography matrix
    Matrix<float> Hinv = inverse(H);

    // Call the helper function to blend the two images
    blendImages(I1, I2, I, Hinv, x0, y0);

    // Display the final blended panorama
    display(I, 0, 0);
}

// Main function
int main(int argc, char* argv[]) {
    const char* s1 = argc>1? argv[1]: srcPath("image0006.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("image0007.jpg");

    // Load and display images
    Image<Color> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load the images" << endl;
        return 1;
    }
    Window w1 = openWindow(I1.width(), I1.height(), s1);
    display(I1,0,0);
    Window w2 = openWindow(I2.width(), I2.height(), s2);
    setActiveWindow(w2);
    display(I2,0,0);

    // Get user's clicks in images
    vector<IntPoint2> pts1, pts2;
    getClicks(w1, w2, pts1, pts2);

    vector<IntPoint2>::const_iterator it;
    cout << "pts1="<<endl;
    for(it=pts1.begin(); it != pts1.end(); it++)
        cout << *it << endl;
    cout << "pts2="<<endl;
    for(it=pts2.begin(); it != pts2.end(); it++)
        cout << *it << endl;

    // Compute homography
    Matrix<float> H = getHomography(pts1, pts2);
    cout << "H=" << H/H(2,2);

    // Apply homography
    panorama(I1, I2, H);

    endGraphics();
    return 0;
}
