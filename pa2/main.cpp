// clang-format off
#include <iostream>
#include <opencv2/opencv.hpp>
#include "rasterizer.hpp"
#include "global.hpp"
#include "Triangle.hpp"
using namespace std;
using namespace cv;
using namespace Eigen;
constexpr double MY_PI = 3.1415926;

Matrix4f get_view_matrix(Vector3f eye_pos)
{
    Matrix4f view = Matrix4f::Identity();

    Matrix4f translate;
    translate << 1,0,0,-eye_pos[0],
                 0,1,0,-eye_pos[1],
                 0,0,1,-eye_pos[2],
                 0,0,0,1;

    view = translate*view;

    return view;
}

Matrix4f get_model_matrix(float rotation_angle)
{
    Matrix4f model = Matrix4f::Identity();
    return model;
}

Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    // Students will implement this function

    Matrix4f projection = Matrix4f::Identity();
    Matrix4f p,s,t;
    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.
    
    p<< zNear,0,0,0, 
        0,zNear,0,0, 
        0,0,zNear+zFar,-zNear*zFar, 
        0,0,1,0;

    float top=zNear*tan(eye_fov/2/180*MY_PI);
    top=-top;   //负号是纠正物体在摄像头后
    float bottom=-top;
    float right=top*aspect_ratio;
    float left=-right;
    
    t<< 1,0,0,-(right+left)/2, 
        0,1,0,-(top+bottom)/2, 
        0,0,1,-(zNear+zFar)/2, 
        0,0,0,1;

    s<< 2/(right-left),0,0,0, 
        0,2/(top-bottom),0,0, 
        0,0,2/(zNear-zFar),0, 
        0,0,0,1;
    projection=s*t*p;
    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = true;
    string filename = "output.png";

    if (argc == 2)
    {
        command_line = true;
        filename = string(argv[1]);
    }

    rst::rasterizer r(700, 700);

    Vector3f eye_pos = {0,0,5};


    vector<Vector3f> pos
            {
                    {2, 0, -2},
                    {0, 2, -2},
                    {-2, 0, -2},
                    {3.5, -1, -5},
                    {2.5, 1.5, -5},
                    {-1, 0.5, -5}
            };

    vector<Vector3i> ind
            {
                    {0, 1, 2},
                    {3, 4, 5}
            };

    vector<Vector3f> cols
            {
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0}
            };

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);
    auto col_id = r.load_colors(cols);

    int key = 0;
    int frame_count = 0;

    if (command_line)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);
        Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cvtColor(image, image, COLOR_RGB2BGR);

        imwrite(filename, image);

        return 0;
    }

    while(key != 27)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);

        Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cvtColor(image, image, COLOR_RGB2BGR);
        imshow("image", image);
        key = waitKey(10);

        cout << "frame count: " << frame_count++ << '\n';
    }

    return 0;
}
// clang-format on