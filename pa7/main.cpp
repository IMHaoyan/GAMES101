#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Sphere.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include <chrono>

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{

    // Change the definition here to change resolution
    Scene scene(784, 784);

    Material* red = new Material(DIFFUSE, Vector3f(0.0f));
    red->Kd = Vector3f(0.63f, 0.065f, 0.05f);
    Material* green = new Material(DIFFUSE, Vector3f(0.0f));
    green->Kd = Vector3f(0.14f, 0.45f, 0.091f);
    Material* white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Material* light = new Material(DIFFUSE, (8.0f * Vector3f(0.747f+0.058f, 0.747f+0.258f, 0.747f) + 15.6f * Vector3f(0.740f+0.287f,0.740f+0.160f,0.740f) + 18.4f *Vector3f(0.737f+0.642f,0.737f+0.159f,0.737f)));
    light->Kd = Vector3f(0.65f);

    float sp = 0.3f;
    Material* Micro_white = new Material(Microfacet, Vector3f(0.0f));
    Micro_white->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Micro_white->Ks = Vector3f(sp,sp,sp);

    // Material* microfacet0 = new Material(Microfacet, Vector3f(0.0f));
   
    // microfacet0->Kd = Vector3f(0.95,0.93,0.88);
    // microfacet0->Ks = Vector3f(sp,sp,sp);

    // Material* microfacet = new Material(Microfacet, Vector3f(0.0f));
    // sp = 0.5f;
    // microfacet->Kd = Vector3f(1.0,0.86,0.57);
    // microfacet->Ks = Vector3f(sp,sp,sp);


    MeshTriangle floor("../models/cornellbox/floor.obj", Micro_white);
    MeshTriangle shortbox("../models/cornellbox/shortbox.obj", Micro_white);
    MeshTriangle tallbox("../models/cornellbox/tallbox.obj", Micro_white);
    MeshTriangle left("../models/cornellbox/left.obj", red);
    MeshTriangle right("../models/cornellbox/right.obj", green);
    MeshTriangle light_("../models/cornellbox/light.obj", light);

    // Sphere sphere((Vector3f(395,100,300)),50.0,microfacet);
    // Sphere sphere2((Vector3f(276,100,300)),50.0,microfacet);
    // Sphere sphere1((Vector3f(157,100,300)),50.0,microfacet0);

    scene.Add(&floor);
    scene.Add(&left);
    scene.Add(&right);
    scene.Add(&light_);
    scene.Add(&tallbox);
    scene.Add(&shortbox);
    // scene.Add(&sphere);
    // scene.Add(&sphere1);
    // scene.Add(&sphere2);
    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();
    r.Render(scene);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}