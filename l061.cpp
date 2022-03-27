// Name: Isabel Wang
// Date: 03/26/2021
// Period: 7
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cmath>
using namespace std;

// global variables
const int dims = 3; // number of dimensions/coordinates/features for each point
int num_points = 50; // number of points in the image, set later in readPPM()
const string input = "image.ppm"; // name of the ppm file
const double high_threshold = 17.0, low_threshold = 11.0, center_threshold = 125.0; // thresholds for strong and weak edges
const double strong = 2.0, weak = 1.0, none = 0.0; // values for edges
const string unvisited = ".", visited = "*"; // markers for unvisited and visited cells
int width = 0, height = 0, maximum = 0, blurred_max = 0; // features of the PPMs, set later in readPPM()

class Point
{
    private:
        int x, y; // x (col) and y (row) coordinates in the image/result matrix
        vector<double> features; // grayscale value
    public:
        Point()
        {
            x = 0;
            y = 0;
            for(int i = 0; i < dims; i++) {
                features.push_back( (double)rand() / RAND_MAX );
            }
        }
        Point(int col, int row, vector<double> coords)
        {
            x = col;
            y = row;
            features = coords;
        }
        ~Point() {}
        double& getFeature(int i) { return features[i]; }
        int getRoundedFeature(int i) { return round(features[i]*num_points); }
        vector<double> getFeatures() { return features; }
        void setFeature(int i, double val) { features[i] = val; }
        int& getX() { return x; }
        int& getY() { return y; }
        double distance(Point&); // finds distance between this point and the given point
        string toString();
};

double Point::distance(Point &p) {
    double total = 0.0;
    for(int i = 0; i < dims; i++) {
        double dcoord = pow(p.getFeature(i) - getFeature(i), 2.0);
        total += dcoord;
    }
    return sqrt(total);
}

string Point::toString() {
    string s;
    int count = 0;
    for(int i = 0; i < features.size(); i++) {
        double x = features[i];
        x = x + 0.5 - (x<0);
        int y = (int)x;
        while(count < dims) {
            s += to_string(y) + " ";
            count++;
        }
        count = 0;
    }
    return s;
}

class Line
{
    private:
        Point p; // line goes through this point
        Point p1, p2; // endpoints of the line
        double m, b; // slope and y-intercept in terms of unit square values
    public:
        Line() {}
        Line(Point &point, double &angle)
        {
            p = point;
            m = tan(angle * M_PI / 180);
            b = (-1 * m * p.getX()) + p.getY(); // y-y0 = m(x-x0) so b = -m*x0 + y0
        }
        ~Line() {}
        double& getSlope() { return m; }
        double& getYIntercept() { return b; }
        Point& getThroughPoint() { return p; }
        Point& getPointOne() { return p1; }
        Point& getPointTwo() { return p2; }
        void extendLine(void);
        void bresenhamsX(vector<Point>&, vector<Point>&, int, int, int, int);
        void bresenhamsY(vector<Point>&, vector<Point>&, int, int, int, int);
        void drawLine(vector<Point>&, vector<Point>&);
        string toString() { return "y=" + to_string(m) + "*x+" + to_string(b); }
};

void Line::bresenhamsX(vector<Point> &pts, vector<Point> &votes, int dx, int dy, int stepx, int stepy) { // x axis is driving axis
    int j = p1.getY(); // j = y1
    int error = 0;
    // find initial error
    if (dx > 0) { error -= dx; }
    else if (dx < 0) { error += dx; }
    if (dy > 0) { error += dy; }
    else if (dy < 0) { error -= dy; }
    // run integer bresenham's algorithm
    for (int i = p1.getX(); (i <= p2.getX() && dx > 0) || (i >= p2.getX() && dx < 0); i += stepx)
    {
        if(i >= 0 && i < width && j >= 0 && j < height)
        {
            int index = j * width + i;
            votes[index].setFeature(0, votes[index].getFeature(0) + 1);
        }
        if (error >= 0)
        {
            j += stepy;
            error -= abs(dx);
        }
        error += abs(dy);
    }
}

void Line::bresenhamsY(vector<Point> &pts, vector<Point> &votes, int dx, int dy, int stepx, int stepy) { // y axis is driving axis
    int i = p1.getX();
    int error = 0;
    // find initial error
    if (dx > 0) { error += dx; }
    else if (dx < 0) { error -= dx; }
    if (dy > 0) { error -= dy; }
    else if (dy < 0) { error += dy; }
    // run integer bresenham's algorithm
    for (int j = p1.getY(); (j <= p2.getY() && dy > 0) || (j >= p2.getY() && dy < 0); j += stepy)
    {
        if(i >= 0 && i < width && j >= 0 && j < height)
        {
            int index = j * width + i;
            votes[index].setFeature(0, votes[index].getFeature(0) + 1);
        }
        if (error >= 0)
        {
            i += stepx;
            error -= abs(dy);
        }
        error += abs(dx);
    }
}

void Line::drawLine(vector<Point> &pts, vector<Point> &votes) {
    // calculate delta x and delta y
    int delta_x = p2.getX() - p1.getX(), delta_y = p2.getY() - p1.getY();
    // determine whether to increment or decrement in x and y directions
    int step_x = 0, step_y = 0;
    if (delta_x > 0) { step_x = 1; }
    else if (delta_x < 0) { step_x = -1; }
    if (delta_y > 0) { step_y = 1; }
    else if (delta_y < 0) { step_y = -1; }
    // determine driving axis and call appropriate method
    if (abs(delta_x) >= abs(delta_y)) { bresenhamsX(pts, votes, delta_x, delta_y, step_x, step_y); }
    else { bresenhamsY(pts, votes, delta_x, delta_y, step_x, step_y); }
}

void Line::extendLine(void) {
    vector<double> v;
    v.push_back(0.0);
    p1 = Point(b, 0, v);
    p2 = Point(m*(width-1)+b, width-1, v);
}

class Circle
{
    public:
        Point c; // center of the circle
        double r; // radius of the circle
    Circle() {}
    Circle(Point center, double radius)
    {
        c = center;
        r = radius;
    }
    void setPixel(int, int, vector<Point>&); // only sets the value in the pixels array for circles, not lines
    void drawCircle(vector<Point>&);
};

void Circle::setPixel(int x, int y, vector<Point> &v) {
    int i = c.getX() + x, j = c.getY() + y; // adjust to circle center
    if(i >= 0 && i < height && j >= 0 && j < width)
    {
        int index = i*width+j;
        v[index].setFeature(0, maximum);
        v[index].setFeature(1, 0.0);
        v[index].setFeature(2, 0.0);
    }
}

void Circle::drawCircle(vector<Point> &v) {
    int x, y, x_max, y2, y2_new, ty;
    x_max = (int)(r / sqrt(2)); // maximum x at radius/sqrt(2)
    y = r;
    y2 = y * y;
    ty = (2 * y) - 1;
    y2_new = y2;
    for (x = 0; x <= x_max; x++)
    {
        if ((y2 - y2_new) >= ty)
        {
            y2 -= ty;
            y -= 1;
            ty -= 2;
        }
        setPixel(x, y, v);
        setPixel(x, -y, v);
        setPixel(-x, y, v);
        setPixel(-x, -y, v);
        setPixel(y, x, v);
        setPixel(y, -x, v);
        setPixel(-y, x, v);
        setPixel(-y, -x, v);
        y2_new -= (2 * x) - 3;
    }
}

void readPPM(vector<Point> &pts) {
    ifstream file;
    file.open(input);
    string in, ppm_type;
    int count = 0; // count = row number in ppm file
    vector<double> values;
    if(file.is_open()) {
        while(getline(file, in)) {
            if(count == 0) {
                ppm_type = in;
            }
            else if(count == 1) {
                width = stoi(in.substr(0, in.find(" ")));
                height = stoi(in.substr(in.find(" ") + 1, in.length() - in.find("  ") - 1));
                num_points = width * height;
            }
            else if(count == 2) {
                maximum = stoi(in);
            }
            else {
                istringstream iss(in); 
                for(string s; iss >> s;)
                    values.push_back(stod(s));
            }
            count++;
        }
    }
    for(int i = 0; i < height; i++) { // iterates through the rows (y coordinate)
        for(int j = 0; j < width; j++) { // iterates through the columns (x coordinate)
            vector<double> c;
            int base = i*width*dims + j*dims;
            for(int l = 0; l < dims; l++) {
                c.push_back(values[base+l]);
            }
            Point p = Point(j, i, c);
            pts.push_back(p);
        }
    }
    file.close();
}

void convertToGrayscale(vector<Point> &pts, vector<Point> &grayscale) {
    // convert rgb values to grayscale value (by averaging them)
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            double total = 0.0;
            for(int l = 0; l < dims; l++) {
                total += pts[i*width + j].getFeature(l);
            }
            double avg_val = total / dims;
            vector<double> c;
            c.push_back(avg_val);
            Point g = Point(j, i, c);
            grayscale.push_back(g);
        }
    }
}

void writePPM(vector<Point> &pts, string filename, int max_val) {
    ofstream file;
    file.open(filename);
    file << "P3" << "\n";
    file << width << " " << height << "\n";
    file << max_val << "\n";
    for (int j = 0; j < height; j++) // keeps track of row/y
    {
        for (int i = 0; i < width; i++) // keeps track of col/x
        {
            if(filename == "imagev2.ppm" && pts[j*width+i].getFeature(0) > center_threshold) {
                file << pts[j*width+i].getFeature(0) << " 0 0 ";
            }
            else if(filename == "imageCC.ppm") {
                file << pts[j*width+i].getFeature(0) << " " << pts[j*width+i].getFeature(1) << " " << pts[j*width+i].getFeature(2) << " ";
            }
            else {
                file << pts[j*width+i].toString();
            }
        }
        file << "\n";
    }
    file.close();
}

void createKernels(vector< vector<int> > &kx, vector< vector<int> > &ky) {
    // vector< vector<int> > kx{{1, 0, -1}, {2, 0, -2}, {1, 0, -1}}, ky{{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}}
    vector<int> vx1, vx2, vx3;
    vx1.push_back(1);
    vx1.push_back(0);
    vx1.push_back(-1);
    vx2.push_back(2);
    vx2.push_back(0);
    vx2.push_back(-2);
    vx3.push_back(1);
    vx3.push_back(0);
    vx3.push_back(-1);
    kx.push_back(vx1);
    kx.push_back(vx2);
    kx.push_back(vx3);
    vector<int> vy1, vy2, vy3;
    vy1.push_back(1);
    vy1.push_back(2);
    vy1.push_back(1);
    vy2.push_back(0);
    vy2.push_back(0);
    vy2.push_back(0);
    vy3.push_back(-1);
    vy3.push_back(-2);
    vy3.push_back(-1);
    ky.push_back(vy1);
    ky.push_back(vy2);
    ky.push_back(vy3);
}

vector<double> getBlurredValue(vector<Point> &grayscale, int index, vector< vector<int> > &kx, vector< vector<int> > &ky) {
    double totalx = 0.0, totaly = 0.0;
    vector<double> gvals;
    for(int j = -1; j < 2; j++) {
        for(int i = -1; i < 2; i++) {
            int curr_index = index + (j*width) + i;
            totalx += (kx[i+1][j+1] * grayscale[curr_index].getFeature(0));
            totaly += (ky[i+1][j+1] * grayscale[curr_index].getFeature(0));
        }
    }
    gvals.push_back(totalx / 9.0);
    gvals.push_back(totaly / 9.0);
    return gvals;
}

void applyOperator(vector<Point> &grayscale, vector<Point> &blurred, vector<double> &thetas, vector< vector<int> > &kx, vector< vector<int> > &ky) {
    vector<double> v;
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            int index = j*width + i;
            if(j == 0 || i == 0 || j == height - 1 || i == width - 1) { // set edges to 0
                v.push_back(0.0);
                blurred.push_back(Point(j, i, v));
                thetas.push_back(0.0);
            }
            else {
                vector<double> bval = getBlurredValue(grayscale, index, kx, ky);
                double g = sqrt(pow(bval[0], 2.0) + pow(bval[1], 2.0));
                double theta = atan2(bval[1], bval[0]) * 180 / M_PI; // angle in degrees (converted from radians)
                v.push_back(g);
                blurred.push_back(Point(j, i, v));
                thetas.push_back(theta);
            }
            v.clear();
        }
    }
}

void transformThetas(vector<double> &thetas, vector<double> &angles, vector<double> &std_angles) {
    for(double x = -180.0; x <= 180.0; x += 45.0) {
        std_angles.push_back(x);
    }
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            double min_diff = 360.0;
            double min_index = 0;
            for(int a = 0; a < std_angles.size(); a++) {
                double diff = abs(std_angles[a] - thetas[j*width+i]);
                if(diff < min_diff || (diff == min_diff && abs(a) > abs(min_index))) {
                    min_diff = diff;
                    min_index = a;
                }
            }
            angles.push_back(std_angles[min_index]);
        }
    }
}

void applyThreshold(vector<Point> &blurred, vector<Point> &edges) {
    vector<double> v;
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            if(blurred[j*width+i].getFeature(0) > high_threshold) {
                v.push_back(strong);
            }
            else if(blurred[j*width+i].getFeature(0) > low_threshold) {
                v.push_back(weak);
            }
            else {
                v.push_back(none);
            }
            edges.push_back(Point(j, i, v));
            v.clear();
        }
    }
}

void areaFill(vector<Point> &edges, vector<Point> &h, vector<string> &record, int r, int c) {
    // base cases: out of bounds, already traversed (mark it as already traversed somehow?), not an edge
    int index = r*width+c;
    if(r >= 0 && r < width && c >= 0 && c < height && (edges[index].getFeature(0) == strong || edges[index].getFeature(0) == weak) && record[index] == unvisited) {
        h[index].setFeature(0, strong);
        record[index] = visited;
        areaFill(edges, h, record, r-1, c-1);
        areaFill(edges, h, record, r-1, c);
        areaFill(edges, h, record, r-1, c+1);
        areaFill(edges, h, record, r, c-1);
        areaFill(edges, h, record, r, c+1);
        areaFill(edges, h, record, r+1, c-1);
        areaFill(edges, h, record, r+1, c);
        areaFill(edges, h, record, r+1, c+1);
    }
}

void applyHysteresis(vector<Point> &edges, vector<Point> &h) {
    vector<string> record; // use this vector to keep track of which spaces have been visited
    // copy edges into h
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            Point p = Point(j, i, edges[j*width+i].getFeatures());
            h.push_back(p);
            record.push_back(unvisited);
        }
    }
    // do area fill when you come across a strong edge
    for(int j = 0; j < height; j++) { // row
        for(int i = 0; i < width; i++) { // column
            if(edges[j*width+i].getFeature(0) == strong) {
                areaFill(edges, h, record, j, i);
            }
        }
    }
    // if a weak edge is not connected to a strong edge (not area filled), then make it a non-edge pixel
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            if(h[j*width+i].getFeature(0) == weak) {
                h[j*width+i].setFeature(0, none);
            }
        }
    }
    // convert to binary (0 or 1), non-edges are already 0s so set strong edges to 1s
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            if(h[j*width+i].getFeature(0) == strong) {
                h[j*width+i].setFeature(0, 1.0);
            }
        }
    }
}

void applyNMS(vector<Point> &edges, vector<Point> &nms, vector<double> &angles) {
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            int index = j * width + i;
            double dx = cos(angles[index] * M_PI / 180), dy = sin(angles[index] * M_PI / 180);
            if(dx != 0 && dx != 1 && dy != 0 && dy != 1) {
               dx *= 2.0 / sqrt(2);
               dy *= 2.0 / sqrt(2);
               dx = round(dx);
               dy = round(dy);
            }
            int i1 = (int) (j + dy) * width + (i + dx);
            int i2 = (int) (j - dy) * width + (i - dx);
            vector<double> v;
            if(i1 >= 0 && i1 < edges.size() && i2 >= 0 && i2 < edges.size()) {
                if(edges[index].getFeature(0) >= edges[i1].getFeature(0) && edges[index].getFeature(0) >= edges[i2].getFeature(0)) {
                    v.push_back(1.0);
                }
                else {
                    v.push_back(0.0);
                }
            }
            else if(i1 >= 0 && i1 < edges.size()) {
                if(edges[index].getFeature(0) >= edges[i1].getFeature(0)) {
                    v.push_back(1.0);
                }
                else {
                    v.push_back(0.0);
                }
            }
            else if(i2 >= 0 && i2 < edges.size()) {
                if(edges[index].getFeature(0) >= edges[i2].getFeature(0)) {
                    v.push_back(1.0);
                }
                else {
                    v.push_back(0.0);
                }
            }
            Point p = Point(j, i, v);
            nms.push_back(p);
            v.clear();
        }
    }
}

void combineNMSH(vector<Point> &final, vector<Point> &nms, vector<Point> &h) {
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            int index = j * width + i;
            vector<double> v;
            if(nms[index].getFeature(0) == weak && h[index].getFeature(0) == weak) {
                v.push_back(1.0);
            }
            else {
                v.push_back(0.0);
            }
            Point p = Point(j, i, v);
            final.push_back(p);
            v.clear();
        }
    }
}

void voting(vector<Point> &votes, vector<Point> &edges, vector<double> &thetas) {
    vector<double> v;
    v.push_back(0);
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            votes.push_back(Point(j, i, v));
        }
    }
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            if(edges[j*width+i].getFeature(0) == 1.0) {
                Line l = Line(edges[j*width+i], thetas[j*width+i]);
                l.extendLine();
                l.drawLine(edges, votes);
            }
        }
    }
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            int count = votes[j*width+i].getFeature(0);
            if(count > blurred_max) {
                blurred_max = count;
            }
        }
    }
}

void findPossibleCenters(vector<Point> &c, vector<Point> &votes, vector<Point> &pts) {
    vector<double> v;
    v.push_back(maximum);
    v.push_back(0);
    v.push_back(0);
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            c.push_back(Point(j, i, pts[j*width+i].getFeatures()));
        }
    }
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            int totalvotes = 0;
            for(int dj = 0; dj < 2; dj++) {
                for(int di = 0; di < 2; di++) {
                    int index = (j+dj)*width+(i+di);
                    if(index < width*height) {
                        totalvotes += votes[index].getFeature(0);
                    }
                }
            }
            if(totalvotes > center_threshold) {
                Point p = Point(j, i, v);
                Circle c1 = Circle(p, 1.0);
                Circle c2 = Circle(p, 2.0);
                Circle c3 = Circle(p, 3.0);
                Circle c4 = Circle(p, 4.0);
                Circle c5 = Circle(p, 5.0);
                c1.drawCircle(c);
                c2.drawCircle(c);
                c3.drawCircle(c);
                c4.drawCircle(c);
                c5.drawCircle(c);
            }
        }
    }
}

void part1() {
    vector<Point> pts; // in pts, the points have a number dims of features
    vector<Point> grayscale_pts; // one feature. this is the image after the original image is converted to grayscale.
    vector<Point> blurred_pts; // one feature. this is the image after the sobel operator is applied.
    vector<double> thetas, angles, std_angles; // thetas for the angles corresponding to G (blurred_pts), angles for "transformed" theta values (in degrees, not radians)
    vector<Point> possibleedges_pts; // one feature. this is the image after the double threshold is applied.
    vector<Point> nms_pts, h_pts; // one feature. these are the images after hysteresis or non-maximum suppression is applied, for canny edge detection
    vector<Point> finaledges_pts; // combination of hysteresis and non-maximum suppression; final, thinned-out edges
    vector<Point> votes, centers;
    readPPM(pts);
    convertToGrayscale(pts, grayscale_pts);
    // writePPM(grayscale_pts, "imageg.ppm", maximum);
    vector< vector<int> > kx, ky; // kernel matrices for the sobel operator, kx will help detect vertical lines, ky will help detect horizontal lines
    createKernels(kx, ky); // add values to the empty kernel matrices
    applyOperator(grayscale_pts, blurred_pts, thetas, kx, ky); // apply the sobel operator
    // writePPM(blurred_pts, "imagem.ppm", blurred_max);
    transformThetas(thetas, angles, std_angles);
    applyThreshold(blurred_pts, possibleedges_pts); // apply double threshold (for hysteresis only, not NMS)
    applyNMS(blurred_pts, nms_pts, angles); // apply non-maximum suppression
    // writePPM(nms_pts, "imagen.ppm", 1);
    applyHysteresis(possibleedges_pts, h_pts); // apply hysteresis (separate from NMS)
    // writePPM(h_pts, "imageh.ppm", 1);
    combineNMSH(finaledges_pts, nms_pts, h_pts); // combine hysteresis & double threshold with non-maximum suppression
    // writePPM(finaledges_pts, "imagef.ppm", 1);
    voting(votes, finaledges_pts, thetas); // result of the voting system
    writePPM(votes, "imagev.ppm", blurred_max);
    // writePPM(votes, "imagev2.ppm", blurred_max);
    findPossibleCenters(centers, votes, pts);
    writePPM(centers, "imageCC.ppm", maximum);
}

int main() {
    part1();
}