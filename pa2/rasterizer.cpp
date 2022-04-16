// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>
using namespace std;
using namespace Eigen;
rst::pos_buf_id rst::rasterizer::load_positions(const vector<Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const vector<Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const vector<Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

static bool insideTriangle(float x, float y, const Vector3f* _v)
{   
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    Vector3f p0p1(_v[0].x() - _v[1].x(), _v[0].y() - _v[1].y(),1.0f);
    Vector3f p1p2(_v[1].x() - _v[2].x(), _v[1].y() - _v[2].y(), 1.0f);
    Vector3f p2p0(_v[2].x() - _v[0].x(), _v[2].y() - _v[0].y(), 1.0f);

    Vector3f p0p(_v[0].x() - x, _v[0].y() - y, 1.0f);
    Vector3f p1p(_v[1].x() - x, _v[1].y() - y, 1.0f);
    Vector3f p2p(_v[2].x() - x, _v[2].y() - y, 1.0f);

    if (p0p1.cross(p0p).z()>0 && p1p2.cross(p1p).z()>0 && p2p0.cross(p2p).z()>0||
    p0p1.cross(p0p).z()<0 && p1p2.cross(p1p).z()<0 && p2p0.cross(p2p).z()<0)
        return true;
    else
        return false;
}
static tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
    float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
    float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
    return {c1,c2,c3};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //Homogeneous division
        for (auto& vec : v) {
            vec /= vec.w();
        }
        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);
    }
    // downsampling 用于SSAA
    /*
    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {
        //每个像素被分成4个小像素，我们在每个小像素上计算和存储depth_buffer和frame_buffer，最后再根据4个小像素的平均值，
        //算出原像素的值。
        //e.g. 假设三角形颜色是红色(255, 0, 0)，有3个小像素在三角形内，那么这个像素的sample list
        //就是(255, 0, 0)，(255, 0, 0)，(255, 0, 0)，(0, 0, 0)，我们通过求这个sample list的平均值得到原像素的值。
           Vector3f color = {0, 0, 0};
            float infinity = numeric_limits<float>::infinity();
            for(float j = start_point; j < 1; j+=pixel_size_sm) {
                for(float i = start_point; i < 1; i+=pixel_size_sm) {
                    int index = get_index_ssaa(x, y, i, j);
                    color += frame_buf_ssaa[index];
                }
            }
            Vector3f p;
            p << x, y, 0.;
            set_pixel(p, color/(ssaa_h*ssaa_w));
        }
    }
    */
}
//SMAA
/*
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    auto v = t.toVector4();

    vector<float> x_arry{ v[0].x(), v[1].x(), v[2].x() };
    vector<float> y_arry{ v[0].y(), v[1].y(), v[2].y() };
    sort(x_arry.begin(), x_arry.end());
    sort(y_arry.begin(), y_arry.end());
    int min_x = floor(x_arry[0]), max_x =ceil( x_arry[2]),
        min_y=floor(y_arry[0]), max_y = ceil(y_arry[2]);
    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle
     for(int y = min_y; y < max_y; y++) {
        for(int x = min_x; x < max_x; x++) {
            // iterate through the sampling point
            for(float j = start_point; j < 1.0; j+=pixel_size_sm) {
                for(float i = start_point; i < 1.0; i+=pixel_size_sm) {
                    // find if the current pixel is inside the triangle
                    if(insideTriangle(x+i, y+j, t.v)) {
                        // if so, use the following code to get the interpolated z value.
                        auto[alpha, beta, gamma] = computeBarycentric2D(x+i, y+j, t.v);
                        float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                        float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                        z_interpolated *= w_reciprocal;

                        // set the current pixel (use the set_pixel function) 
                        // to the color of the triangle (use getColor function)
                        // if it should be painted.
                        int index = get_index_ssaa(x, y, i, j);
                        if(z_interpolated < depth_buf_ssaa[index]) {
                            frame_buf_ssaa[index] = t.getColor();
                            depth_buf_ssaa[index] = z_interpolated;
                        }
                    }
                }
            }
        }
    }
}
*/

//无抗锯齿
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    auto v = t.toVector4();
    float min_x = width;
    float max_x = 0;
    float min_y = height;
    float max_y = 0;
    // find out the bounding box of current triangle
    for(int i = 0; i < 3; i++) {
        min_x = std::min(v[i].x(), min_x);
        max_x = std::max(v[i].x(), max_x);
        min_y = std::min(v[i].y(), min_y);
        max_y = std::max(v[i].y(), max_y);
    }
    // iterate through the pixel
    for(int y = min_y; y < max_y; y++) {
        for(int x = min_x; x < max_x; x++) {
            int index = get_index(x, y);
            // find if the current pixel is inside the triangle
            if(insideTriangle(x+0.5, y+0.5, t.v)) {
                // if so, use the following code to get the interpolated z value.
                auto[alpha, beta, gamma] = computeBarycentric2D(x+0.5, y+0.5, t.v);
                float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                z_interpolated *= w_reciprocal;
                // set the current pixel (use the set_pixel function) 
                // to the color of the triangle (use getColor function)
                // if it should be painted.
                if(z_interpolated < depth_buf[index]) {
                    Vector3f p;
                    p << x, y, z_interpolated;
                    Vector3f color = t.getColor();
                    set_pixel(p, color);
                    depth_buf[index] = z_interpolated;
                }   
            }
        }
    }
}
//MSAA对于每个像素，我们通过三角形在这个像素的覆盖率，即它的面积，来计算原像素的值。我们先把这个像素分成4个小像素，然后
//用有多少个小像素在三角形内，来近似这个面积。e.g. 假设三角形颜色是红色(255, 0, 0)，有3个小像素在三角形内，那么这个像素
//最终的值就是3/4乘以(255, 0, 0)。
/*
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    auto v = t.toVector4();

    float min_x = width;
    float max_x = 0;
    float min_y = height;
    float max_y = 0;

    // find out the bounding box of current triangle
    for(int i = 0; i < 3; i++) {
        min_x = std::min(v[i].x(), min_x);
        max_x = std::max(v[i].x(), max_x);
        min_y = std::min(v[i].y(), min_y);
        max_y = std::max(v[i].y(), max_y);
    }

    // iterate through the pixel
    for(int y = min_y; y < max_y; y++) {
        for(int x = min_x; x < max_x; x++) {
            int index = get_index(x, y);
            float count = 0.0;
            float max_count = ssaa_w*ssaa_h;
            // iterate through the sampling points
            for(float j = start_point; j < 1.0; j+=pixel_size_sm) {
                for(float i = start_point; i < 1.0; i+=pixel_size_sm) {
                    if(insideTriangle(x+i, y+j, t.v)) {
                        count += 1.0;
                    }
                }
            }
            // find if the current pixel is inside the triangle
            if(insideTriangle(x+0.5, y+0.5, t.v)) {
                // if so, use the following code to get the interpolated z value.
                auto[alpha, beta, gamma] = computeBarycentric2D(x+0.5, y+0.5, t.v);
                float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                z_interpolated *= w_reciprocal;

                // set the current pixel (use the set_pixel function) 
                // to the color of the triangle (use getColor function)
                // if it should be painted.
                if(z_interpolated < depth_buf[index]) {
                    Vector3f p;
                    p << x, y, z_interpolated;
                    Vector3f color = t.getColor()*(count/max_count);
                    set_pixel(p, color);
                    depth_buf[index] = z_interpolated;
                }   
            }
        }
    }
}*/
void rst::rasterizer::set_model(const Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        fill(frame_buf.begin(), frame_buf.end(), Vector3f{0, 0, 0});
        fill(frame_buf_ssaa.begin(), frame_buf_ssaa.end(), Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        fill(depth_buf.begin(), depth_buf.end(), numeric_limits<float>::infinity());
        fill(depth_buf_ssaa.begin(), depth_buf_ssaa.end(), numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
    frame_buf_ssaa.resize(w * ssaa_w * h * ssaa_h);
    depth_buf_ssaa.resize(w * ssaa_w * h * ssaa_h);
}
// for float args
int rst::rasterizer::get_index_ssaa(int x, int y, float i, float j)
{
    int ssaa_height = height * ssaa_h;
    int ssaa_width = width * ssaa_w;

    i = int((i-start_point)/pixel_size_sm);
    j = int((j-start_point)/pixel_size_sm);

    return (ssaa_height-1-y*ssaa_h+j)*ssaa_width + x*ssaa_w + i;
}
///
int rst::rasterizer::get_index(int x, int y)
{
    return (height-1-y)*width + x;
}

void rst::rasterizer::set_pixel(const Vector3f& point, const Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height-1-point.y())*width + point.x();
    frame_buf[ind] = color;

}

// clang-format on