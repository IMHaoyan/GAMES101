#include <chrono>
#include <iostream>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;
vector<Point2f> control_points;

void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == EVENT_LBUTTONDOWN && control_points.size() < 4) 
    {
        cout << "Left button of the mouse is clicked - position (" << x << ", "
        << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}

void naive_bezier(const vector<Point2f> &points, Mat &window) 
{
    auto &p_0 = points[0];
    auto &p_1 = points[1];
    auto &p_2 = points[2];
    auto &p_3 = points[3];

    for (double t = 0.0; t <= 1.0; t += 0.001) 
    {
        auto point = pow(1 - t, 3) * p_0 + 3 * t * pow(1 - t, 2) * p_1 +
                 3 * pow(t, 2) * (1 - t) * p_2 + pow(t, 3) * p_3;

        window.at<Vec3b>(point.y, point.x)[2] = 255;
    }
}

Point2f recursive_bezier(const vector<Point2f> &control_points, float t) 
{
    // TODO: Implement de Casteljau's algorithm
    if(control_points.size()==1)return control_points[0];
    vector<Point2f> next_control_points = {};
    for(int i = 0; i < control_points.size()-1; i++) {
        auto &a = control_points[i];
        auto &b = control_points[i+1];
        auto p = (1-t)*a + t*b;
        next_control_points.push_back(p);
    }
    return recursive_bezier(next_control_points,t);
}

void bezier(const vector<Point2f> &control_points, Mat &window) 
{
    // TODO: Iterate through all t = 0 to t = 1 with small steps, and call de Casteljau's 
    // recursive Bezier algorithm.
    for (double t = 0.0; t <= 1.0; t += 0.001) 
    {
        auto p = recursive_bezier(control_points,t);
        //window.at<Vec3b>(point.y, point.x)[1] = 255;
        Point2f p0(p.x-floor(p.x) < 0.5 ? floor(p.x) : ceil(p.x),
                    p.y-floor(p.y) < 0.5 ? floor(p.y) : ceil(p.y));
        vector<Point2f> ps = { p0, Point2f(p0.x-1, p0.y),
            Point2f(p0.x, p0.y-1), Point2f(p0.x-1, p0.y-1),
        };

        float sum_d = 0.f;
        float max_d = sqrt(2);
        vector<float> ds = {};
        for (int i = 0; i < 4; i++) {
            Point2f cp(ps[i].x + 0.5f, ps[i].y + 0.5f);
            float d = max_d - sqrt(pow(p.x - cp.x, 2) + pow(p.y - cp.y, 2));
            ds.push_back(d);
            sum_d += d;
        };

        for (int i = 0; i < 4; i++) {
            float k = ds[i]/sum_d;
            window.at<Vec3b>(ps[i].y, ps[i].x)[1] = min(255.f, window.at<Vec3b>(ps[i].y, ps[i].x)[1] + 255.f * k);
        };
    }
}

int main() 
{
    Mat window = Mat(700, 700, CV_8UC3, Scalar(0));
    cvtColor(window, window, COLOR_BGR2RGB);
    namedWindow("Bezier Curve", WINDOW_AUTOSIZE);

    setMouseCallback("Bezier Curve", mouse_handler, nullptr);

    int key = -1;
    while (key != 27) 
    {
        for (auto &point : control_points) 
        {
            circle(window, point, 3, {255, 255, 255}, 3);
        }

        if (control_points.size() == 4) 
        {
            naive_bezier(control_points, window);
            bezier(control_points, window);

            imshow("Bezier Curve", window);
            imwrite("my_bezier_curve.png", window);
            key = waitKey(0);

            return 0;
        }

        imshow("Bezier Curve", window);
        key = waitKey(20);
    }

    return 0;
}
