//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"

enum MaterialType { DIFFUSE, SPECULAR, GLASS, Microfacet };

class Material {
   private:
    // Compute reflection direction
    Vector3f reflect(const Vector3f& I, const Vector3f& N) const {
        return I - 2 * dotProduct(I, N) * N;
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and
    // negate the normal N
    Vector3f refract(const Vector3f& I,
                     const Vector3f& N,
                     const float& ior) const {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) {
            cosi = -cosi;
        } else {
            swap(etai, etat);
            n = -N;
        }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    
    Vector3f toWorld(const Vector3f& a, const Vector3f& N) {
        Vector3f B, C;
        if (fabs(N.x) > fabs(N.y)) {
            float invLen = 1.0f / sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x * invLen);
        } else {
            float invLen = 1.0f / sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y * invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

   public:
    MaterialType m_type;
    // Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd, Ks;
    float metallic = 0.0;
    float specularExponent;
    // Texture tex;

    inline Material(MaterialType t = DIFFUSE, Vector3f e = Vector3f(0, 0, 0));
    inline MaterialType getType();
    // inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f& wi, const Vector3f& N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f& wi, const Vector3f& wo, const Vector3f& N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f& wi,
                         const Vector3f& wo,
                         const Vector3f& N);

    float fresnel(const Vector3f& I,
                 const Vector3f& N,
                 const float& ior) const {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {
            swap(etai, etat);
        }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            return 1.0;
        } else {
            float cost = sqrtf(max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) /
                       ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) /
                       ((etai * cosi) + (etat * cost));
            return (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is
        // given by: kt = 1 - kr;
    }

    float DistributionGGX(Vector3f N, Vector3f H, float roughness) {
        float a = roughness * roughness;
        float a2 = a * a;
        float NdotH = max(dotProduct(N, H), 0.0f);
        float NdotH2 = NdotH * NdotH;

        float nom = a2;
        float denom = (NdotH2 * (a2 - 1.0) + 1.0);
        denom = M_PI * denom * denom;

        return nom /max(denom, 0.0000001f);  
        // prevent divide by zero for
        // roughness=0.0 and NdotH=1.0
    }

    float SchlickGGX(float NdotV, float roughness) {
        float r = (roughness + 1.0);
        float k = (r * r) / 8.0;

        float nom = NdotV;
        float denom = NdotV * (1.0 - k) + k;

        return nom / denom;
    }

    float GeometrySmith(Vector3f N, Vector3f V, Vector3f L, float roughness) {
        float NdotV = max(dotProduct(N, V), 0.0f);
        float NdotL = max(dotProduct(N, L), 0.0f);
        float ggx2 = SchlickGGX(NdotV, roughness);
        float ggx1 = SchlickGGX(NdotL, roughness);

        return ggx1 * ggx2;
    }
};

Material::Material(MaterialType t, Vector3f e) {
    m_type = t;
    // m_color = c;
    m_emission = e;
}

MaterialType Material::getType() {
    return m_type;
}
/// Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {
    return m_emission;
}
bool Material::hasEmission() {
    if (m_emission.norm() > EPSILON)
        return true;
    else
        return false;
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}

Vector3f Material::sample(const Vector3f& wi, const Vector3f& N) {
    switch (m_type) {
        case DIFFUSE: {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = fabs(1.0f - 2.0f * x_1);
            float r = sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r * cos(phi), r * sin(phi), z);
            return toWorld(localRay, N);

            break;
        }
        case SPECULAR: {
            // uniform sample on the hemisphere
            return reflect(wi, N);
            break;
        }
        case GLASS: {
            return refract(wi, N, 3);
            break;
        }
        case Microfacet: {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = fabs(1.0f - 2.0f * x_1);
            float r = sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r * cos(phi), r * sin(phi), z);
            return toWorld(localRay, N);

            break;
        }
    }
    return Vector3f();
}

float Material::pdf(const Vector3f& wi, const Vector3f& wo, const Vector3f& N) {
    switch (m_type) {
        case DIFFUSE: {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
        case SPECULAR: {
            return 1.0f;
            break;
        }
        case GLASS: {
            return 1.0f;
            break;
        }
        case Microfacet: {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
    }
    return 0.0f;
}

Vector3f Material::eval(const Vector3f& wi,
                        const Vector3f& wo,
                        const Vector3f& N) {
    switch (m_type) {
        case DIFFUSE: {
            // calculate the contribution of diffuse   model
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            } else
                return Vector3f(0.0f);
            break;
        }
        case SPECULAR: {
            return Vector3f(1.0f, 1.0f, 1.0f);
            break;
        }
        case GLASS: {
            return Vector3f(1.0f, 1.0f, 1.0f);
            break;
        }
        case Microfacet: {
            // Disney PBR 方案
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                float roughness = 0.35, etat =1.85;

                Vector3f V = -wi;
                Vector3f L = wo;
                Vector3f H = normalize(V + L);
                float D = DistributionGGX(N, H, roughness);
                float G = GeometrySmith(N, V, L, roughness);
                float F= fresnel(wi, N, etat);

                Vector3f nominator = D * G * F;
                float denominator = 4 * max(dotProduct(N, V), 0.0f) *
                                    max(dotProduct(N, L), 0.0f);
                Vector3f specular = nominator / max(denominator, 0.001f);

                // 能量守恒
                float ks_ = F;
                float kd_ = 1.0f - ks_;

                Vector3f diffuse = (1.0f-metallic) / M_PI;

                // 因为在 specular
                // 项里已经考虑了折射部分的比例：F，所以折射部分不需要再乘以 ks_
                // （ks_ * Ks * specular）
                //return kd_ * Kd * diffuse + Ks * specular;
                return (1.0 - F) * Kd * diffuse + Ks * specular;
            } else
                return Vector3f(0.0f);
            break;
        }
    }
    return Vector3f();
}

#endif  // RAYTRACING_MATERIAL_H
