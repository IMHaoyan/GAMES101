#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;
using namespace Eigen;
constexpr double MY_PI = 3.1415926;

Matrix4f get_view_matrix(Vector3f eye_pos)
{
    Matrix4f view = Matrix4f::Identity();

    Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

Matrix4f get_model_matrix(float rotation_angle)
{
    Matrix4f model = Matrix4f::Identity();

    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.
    float delta=rotation_angle;
    float fa=acos(-1)*delta/180;
    model<<cos(fa),-sin(fa),0,0,sin(fa),cos(fa),0,0,0,0,1,0,0,0,0,1;
    return model;
}

Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    // Students will implement this function

    Matrix4f projection = Matrix4f::Identity();
    Matrix4f p,o,n;
    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.
    p<< zNear,0,0,0, 
        0,zNear,0,0, 
        0,0,zNear+zFar,-zNear*zFar, 
        0,0,1,0;
    float t=tan(eye_fov/2/180*acos(-1))*zNear;
    float b=-t;
    float r=t*aspect_ratio;
    float l=-r;
    o<< 2/(r-l),0,0,0, 
        0,2/(t-b),0,0, 
        0,0,2/(zNear-zFar),0, 
        0,0,0,1;
    n<< 1,0,0,-(r+l)/2, 
        0,1,0,-(t+b)/2, 
        0,0,1,-(zNear+zFar)/2, 
        0,0,0,1;
    projection=o*n*p;
    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = string(argv[3]);
        }
        else
            return 0;
    }

    rst::rasterizer r(700, 700);

    Vector3f eye_pos = {0, 0, 5};

    vector<Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};

    vector<Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        imshow("image", image);
        key = waitKey(10);

        cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}
