//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"
using namespace std;

void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    Intersection p_inter = intersect(ray);
    if (!p_inter.happened) {
        return Vector3f();
    }
    if (p_inter.m->hasEmission()) {//光源处
        return p_inter.m->getEmission();
    }

    float EPLISON = 0.0001;
    //1.直接光照. 
    Intersection x_inter;float pdf_light;
    //找到(一个)光源
    sampleLight(x_inter, pdf_light);
    // Get x, ws, NN, emit from inter
    Vector3f p = p_inter.coords, x = x_inter.coords;
    Vector3f ws_dir = (x - p).normalized();
    float ws_distance = (x - p).norm();
    Vector3f N = p_inter.normal.normalized(), NN = x_inter.normal.normalized();
    Vector3f emit = x_inter.emit;
    Vector3f l_dir=0.f;
    Ray ws_Ray(p, ws_dir);
    Intersection ws_ray_inter = intersect(ws_Ray);
    if(ws_ray_inter.distance - ws_distance > -EPLISON){//光源未被遮挡
        l_dir = emit * p_inter.m->eval(ray.direction,ws_Ray.direction, N) * dotProduct(N, ws_Ray.direction) 
        * dotProduct(NN, -ws_Ray.direction) / pow(ws_distance, 2) / pdf_light;
        //此时cos(theta')是因为进行了积分区间的变换，由出发点的w到光源的A
    }
    float _RussianRoulette = depth > maxDepth ? RussianRoulette: 1.0;
    //float _RussianRoulette = RussianRoulette;
    //2.间接光照
    Vector3f l_indir=0.f;
    if(get_random_float() > _RussianRoulette) {
        return l_dir;//直接返回 此时 l_indir为 0
    }
    Vector3f wi_dir = p_inter.m->sample(ray.direction, N).normalized();
    Ray wi_Ray(p,wi_dir);
    Intersection wi_ray_inter = intersect(wi_Ray);
    float pdf = p_inter.m->pdf(ray.direction, wi_Ray.direction, N);
    if(pdf < EPLISON){
        return l_dir;
    }
    if(wi_ray_inter.happened && !wi_ray_inter.m->hasEmission()){
        l_indir = castRay(wi_Ray, depth+1) * p_inter.m->eval(ray.direction,wi_Ray.direction, N)
            *dotProduct(N, wi_Ray.direction) 
            / pdf / _RussianRoulette;
    }
    return l_dir + l_indir;
}