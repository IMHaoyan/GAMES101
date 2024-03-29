#include <iostream>
#include <opencv2/opencv.hpp>

#include "global.hpp"
#include "rasterizer.hpp"
#include "Triangle.hpp"
#include "Shader.hpp"
#include "Texture.hpp"
#include "OBJ_Loader.h"
using namespace std;
using namespace Eigen;
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

Matrix4f get_model_matrix(float angle)
{
    Matrix4f rotation;
    angle = angle * MY_PI / 180.f;
    rotation << cos(angle), 0, sin(angle), 0,
                0, 1, 0, 0,
                -sin(angle), 0, cos(angle), 0,
                0, 0, 0, 1;

    Matrix4f scale;
    scale << 2.5, 0, 0, 0,
              0, 2.5, 0, 0,
              0, 0, 2.5, 0,
              0, 0, 0, 1;

    Matrix4f translate;
    translate << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

    return translate * rotation * scale;
}

Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    // TODO: Use the same projection matrix from the previous assignments
    Matrix4f projection = Matrix4f::Identity();
    Matrix4f p,s,t;
    
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

Vector3f vertex_shader(const vertex_shader_payload& payload)
{
    return payload.position;
}

Vector3f normal_fragment_shader(const fragment_shader_payload& payload)
{
    Vector3f return_color = (payload.normal.head<3>().normalized() + Vector3f(1.0f, 1.0f, 1.0f)) / 2.f;
    Vector3f result;
    result << return_color.x() * 255, return_color.y() * 255, return_color.z() * 255;
    return result;
}

static Vector3f reflect(const Vector3f& vec, const Vector3f& axis)
{
    auto costheta = vec.dot(axis);
    return (2 * costheta * axis - vec).normalized();
}

struct light
{
    Vector3f position;
    Vector3f intensity;
};

Vector3f texture_fragment_shader(const fragment_shader_payload& payload)
{
    Vector3f return_color = {0, 0, 0};
    if (payload.texture)
    {
        // TODO: Get the texture value at the texture coordinates of the current fragment
        return_color = payload.texture->getColor(payload.tex_coords.x(),payload.tex_coords.y());
    }
    Vector3f texture_color;
    texture_color << return_color.x(), return_color.y(), return_color.z();

    Vector3f kd = texture_color / 255.f;
    Vector3f ka = Vector3f(0.005, 0.005, 0.005);
    Vector3f ks = Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20, 20, 20}, {500, 500, 500}};
    auto l2 = light{{-20, 20, 0}, {500, 500, 500}};

    vector<light> lights = {l1, l2};
    Vector3f amb_light_intensity{10, 10, 10};
    Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Vector3f color = texture_color;
    Vector3f point = payload.view_pos;
    Vector3f normal = payload.normal;

    Vector3f result_color = {0, 0, 0};

    for (auto& light : lights)
    {
        // TODO: For each light source in the code, calculate what the *ambient*, *diffuse*, and *specular* 
        // components are. Then, accumulate that result on the *result_color* object.
        Vector3f l = (light.position - point).normalized();
        Vector3f v = (eye_pos - point).normalized();
        Vector3f h = (v + l).normalized();
        float r2 = (light.position - point).squaredNorm();

        // NOTE: use cwiseProduct() because ka/kd/ks are coefficients
        Vector3f diffuse = kd.cwiseProduct(light.intensity/r2)*max(0.f, normal.dot(l));
        Vector3f specular = ks.cwiseProduct(light.intensity/r2)*pow(max(0.f, normal.dot(h)), p);
        Vector3f ambient = ka.cwiseProduct(amb_light_intensity);
        result_color += (diffuse + specular + ambient);
    }
    return result_color * 255.f;
}

Vector3f phong_fragment_shader(const fragment_shader_payload& payload)
{
    Vector3f kd = payload.color;
    Vector3f ka = Vector3f(0.005, 0.005, 0.005);
    Vector3f ks = Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20, 20, 20}, {500, 500, 500}};
    auto l2 = light{{-20, 20, 0}, {500, 500, 500}};

    vector<light> lights = {l1, l2};
    Vector3f amb_light_intensity{10, 10, 10};
    Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Vector3f color = payload.color;
    Vector3f point = payload.view_pos;
    Vector3f normal = payload.normal;//normal line:法线

    Vector3f result_color = {0, 0, 0};
    for (auto& light : lights)
    {
        // TODO: For each light source in the code, calculate what the *ambient*, *diffuse*, and *specular* 
        // components are. Then, accumulate that result on the *result_color* object.
        Vector3f l = (light.position - point).normalized();
        Vector3f v = (eye_pos - point).normalized();
        Vector3f h = (v + l).normalized();
        float r2 = (light.position - point).squaredNorm();

        // NOTE: use cwiseProduct() because ka/kd/ks are coefficients
        Vector3f diffuse = kd.cwiseProduct(light.intensity/r2)*max(0.f, normal.dot(l));
        Vector3f specular = ks.cwiseProduct(light.intensity/r2)*pow(max(0.f, normal.dot(h)), p);
        Vector3f ambient = ka.cwiseProduct(amb_light_intensity);
        result_color += (diffuse + specular + ambient);
    }
    return result_color * 255.f;
}

Vector3f displacement_fragment_shader(const fragment_shader_payload& payload)
{
    
    Vector3f ka = Vector3f(0.005, 0.005, 0.005);
    Vector3f kd = payload.color;
    Vector3f ks = Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20, 20, 20}, {500, 500, 500}};
    auto l2 = light{{-20, 20, 0}, {500, 500, 500}};

    vector<light> lights = {l1, l2};
    Vector3f amb_light_intensity{10, 10, 10};
    Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Vector3f color = payload.color; 
    Vector3f point = payload.view_pos;
    Vector3f normal = payload.normal;

    float kh = 0.2, kn = 0.1;
    
    // TODO: Implement displacement mapping here
    Vector3f result_color = {0, 0, 0};
    Vector3f n=normal;
    float x=n.x(),y=n.y(),z=n.z();
    Vector3f t(x*y/sqrt(x*x+z*z),sqrt(x*x+z*z),z*y/sqrt(x*x+z*z));
    Vector3f b = n.cross(t);
    Matrix3f TBN;TBN<<t,b,n;    //注意此处为Matrix！！
    float u = payload.tex_coords.x(), v = payload.tex_coords.y();
    int w = payload.texture->width, h = payload.texture->height;
    auto dU = kh * kn * (payload.texture->getColor(u+1.f/w,v).norm()-payload.texture->getColor(u,v).norm());
    auto dV = kh * kn * (payload.texture->getColor(u,v+1.f/h).norm()-payload.texture->getColor(u,v).norm());
    Vector3f ln(-dU, -dV, 1.f);
    n = (TBN * ln).normalized();
    point = point + kn * n * payload.texture->getColor(u,v).norm();// Position p = p + kn * n * h(u,v)

    for (auto& light : lights){
        // TODO: For each light source in the code, calculate what the *ambient*, *diffuse*, and *specular* 
        // components are. Then, accumulate that result on the *result_color* object.
        Vector3f l = (light.position - point).normalized();
        Vector3f v = (eye_pos - point).normalized();
        Vector3f h = (v + l).normalized();
        float r2 = (light.position - point).squaredNorm();
        Vector3f diffuse = kd.cwiseProduct(light.intensity/r2)*max(0.f, normal.dot(l));
        Vector3f specular = ks.cwiseProduct(light.intensity/r2)*pow(max(0.f, normal.dot(h)), p);
        Vector3f ambient = ka.cwiseProduct(amb_light_intensity);
        result_color += (diffuse + specular + ambient);
    }
    
    return result_color * 255.f;
}

Vector3f bump_fragment_shader(const fragment_shader_payload& payload)
{
    Vector3f kd = payload.color;
    Vector3f ks = Vector3f(0.7937, 0.7937, 0.7937);
    Vector3f ka = Vector3f(0.005, 0.005, 0.005);

    auto l1 = light{{20, 20, 20}, {500, 500, 500}};
    auto l2 = light{{-20, 20, 0}, {500, 500, 500}};

    vector<light> lights = {l1, l2};
    Vector3f amb_light_intensity{10, 10, 10};
    Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Vector3f color = payload.color; 
    Vector3f point = payload.view_pos;
    Vector3f normal = payload.normal;

    float kh = 0.2, kn = 0.1;
    // TODO: Implement bump mapping here
    // Let n = normal = (x, y, z)
    // Vector t = (x*y/sqrt(x*x+z*z),sqrt(x*x+z*z),z*y/sqrt(x*x+z*z))
    // Vector b = n cross product t
    // Matrix TBN = [t b n]
    // dU = kh * kn * (h(u+1/w,v)-h(u,v))
    // dV = kh * kn * (h(u,v+1/h)-h(u,v))
    // Vector ln = (-dU, -dV, 1)
    // Normal n = normalize(TBN * ln)
    Vector3f result_color = {0, 0, 0};

    Vector3f n=normal;
    float x=n.x(),y=n.y(),z=n.z();
    Vector3f t(x*y/sqrt(x*x+z*z),sqrt(x*x+z*z),z*y/sqrt(x*x+z*z));
    Vector3f b = n.cross(t);
    Matrix3f TBN;TBN<<t,b,n;    //注意此处为Matrix！！
    float u = payload.tex_coords.x(), v = payload.tex_coords.y();
    int w = payload.texture->width, h = payload.texture->height;
    auto dU = kh * kn * (payload.texture->getColor(u+1.f/w,v).norm()-payload.texture->getColor(u,v).norm());
    auto dV = kh * kn * (payload.texture->getColor(u,v+1.f/h).norm()-payload.texture->getColor(u,v).norm());
    Vector3f ln(-dU, -dV, 1.f);
    n = (TBN * ln).normalized();
    
    result_color = n;
    return result_color * 255.f;
}

int main(int argc, const char** argv)
{
    vector<Triangle*> TriangleList;

    float angle = 140.0;
    bool command_line = false;

    string filename = "output.png";
    objl::Loader Loader;
    string obj_path = "../models/spot/";

    // Load .obj File
    bool loadout = Loader.LoadFile("../models/spot/spot_triangulated_good.obj");
    for(auto mesh:Loader.LoadedMeshes)
    {
        for(int i=0;i<mesh.Vertices.size();i+=3)
        {
            Triangle* t = new Triangle();
            for(int j=0;j<3;j++)
            {
                t->setVertex(j,Vector4f(mesh.Vertices[i+j].Position.X,mesh.Vertices[i+j].Position.Y,mesh.Vertices[i+j].Position.Z,1.0));
                t->setNormal(j,Vector3f(mesh.Vertices[i+j].Normal.X,mesh.Vertices[i+j].Normal.Y,mesh.Vertices[i+j].Normal.Z));
                t->setTexCoord(j,Vector2f(mesh.Vertices[i+j].TextureCoordinate.X, mesh.Vertices[i+j].TextureCoordinate.Y));
            }
            TriangleList.push_back(t);
        }
    }

    rst::rasterizer r(700, 700);

    auto texture_path = "hmap.jpg";
    r.set_texture(Texture(obj_path + texture_path));

    function<Vector3f(fragment_shader_payload)> active_shader = phong_fragment_shader;

    if (argc >= 2)
    {
        command_line = true;
        filename = string(argv[1]);

        if (argc == 3 && string(argv[2]) == "texture")
        {
            filename = "texture.png";
            cout << "Rasterizing using the texture shader\n";
            active_shader = texture_fragment_shader;
            texture_path = "spot_texture.png";
            r.set_texture(Texture(obj_path + texture_path));
        }
        else if (argc == 3 && string(argv[2]) == "normal")
        {
            filename = "normal.png";
            cout << "Rasterizing using the normal shader\n";
            active_shader = normal_fragment_shader;
        }
        else if (argc == 3 && string(argv[2]) == "phong")
        {
            filename = "phong.png";
            cout << "Rasterizing using the phong shader\n";
            active_shader = phong_fragment_shader;
        }
        else if (argc == 3 && string(argv[2]) == "bump")
        {
            filename = "bump.png";
            cout << "Rasterizing using the bump shader\n";
            active_shader = bump_fragment_shader;
        }
        else if (argc == 3 && string(argv[2]) == "displacement")
        {
            filename = "displacement.png";
            cout << "Rasterizing using the bump shader\n";
            active_shader = displacement_fragment_shader;
        }
    }

    Vector3f eye_pos = {0,0,10};

    r.set_vertex_shader(vertex_shader);
    r.set_fragment_shader(active_shader);

    int key = 0;
    int frame_count = 0;

    if (command_line)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);
        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45.0, 1, 0.1, 50));

        r.draw(TriangleList);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imwrite(filename, image);

        return 0;
    }

    while(key != 27)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45.0, 1, 0.1, 50));

        //r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);
        r.draw(TriangleList);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imshow("image", image);
        cv::imwrite(filename, image);
        key = cv::waitKey(10);

        if (key == 'a' )
        {
            angle -= 0.1;
        }
        else if (key == 'd')
        {
            angle += 0.1;
        }

    }
    return 0;
}
