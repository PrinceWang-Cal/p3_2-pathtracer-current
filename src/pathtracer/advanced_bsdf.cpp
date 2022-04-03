#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Project 3-2: Part 1
  // Implement MirrorBSDF
    reflect(wo, wi);
    *pdf = 1.0;
  return (this->reflectance) / abs_cos_theta(*wi);
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
  // TODO Project 3-2: Part 2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
  Vector3D n = Vector3D(0, 0, 1);
  double cos_theta_h = dot(h, n)/ (h.norm() * n.norm());
  double sin2_theta_h = 1 - pow(cos_theta_h, 2);
  double tan2_theta_h = sin2_theta_h / pow(cos_theta_h, 2);
    
  double numerator = exp(-1 * tan2_theta_h / pow(alpha, 2));
    double denomenator = PI * pow(alpha, 2) * pow(cos_theta_h, 4);
  return numerator/denomenator;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO Project 3-2: Part 2
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.
    Vector3D r_s_numerator = (eta * eta + k * k) - 2.0 * eta * cos_theta(wi) + cos_theta(wi) * cos_theta(wi);
    Vector3D r_s_denomenator = (eta * eta + k * k) + 2.0 * eta * cos_theta(wi) + cos_theta(wi) * cos_theta(wi);
    
    Vector3D r_p_numerator = (eta * eta + k * k) * cos_theta(wi) * cos_theta(wi) -  2.0 * eta * cos_theta(wi) + 1.0;
    Vector3D r_p_denomenator = (eta * eta + k * k) * cos_theta(wi) * cos_theta(wi) +  2.0 * eta * cos_theta(wi) + 1.0;
    
    Vector3D r_s = r_s_numerator / r_s_denomenator;
    Vector3D r_p = r_p_numerator / r_p_denomenator;

    return (r_s + r_p) / 2.0;
}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO Project 3-2: Part 2
  // Implement microfacet model here.
    Vector3D n = Vector3D(0, 0, 1);
    if (dot(n, wo) < 0) {
          return 0;
    } else if (dot(n, wi) < 0) {
          return 0;
    }
    
    Vector3D h = (wo + wi) / (wo + wi).norm();
    Vector3D numerator = F(wi) * G(wo, wi) * D(h);
    Vector3D denomenator = 4 * dot(n, wo) * dot(n, wi);

    return numerator/denomenator;
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Project 3-2: Part 2
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
    
    Vector2D r = sampler.get_sample();
    double r1 = r[0];
    double r2 = r[1];
    Vector3D n = Vector3D(0, 0, 1);

    double theta_h = atan(sqrt(-1 * pow(alpha, 2) * log(1 - r1)));
    double phi_h = 2.0 * PI * r2;
    double h1 = sin(theta_h) * cos(phi_h);
    double h2 = sin(theta_h) * sin(phi_h);
    double h3 = cos(theta_h);

    Vector3D h = Vector3D(h1,h2,h3);

    *wi = -1.0 * wo + 2.0 * dot(h, wo) * h;

    //calculate pdf
    if (dot(n, *wi) < 0) {
          *pdf = 0;
          return 0;
    } else {
        double p_theta_1 = 2 * sin(theta_h);
        double p_theta_2 = exp(-1 * pow(tan(theta_h), 2) / pow(alpha, 2));
        double p_theta_3 = (pow(alpha, 2) * pow(cos(theta_h), 3));


        double p_theta = (p_theta_1 * p_theta_2) / (p_theta_3);
        double p_phi = 1 / (2 * PI);
        double p_h = p_theta * p_phi / sin(theta_h);

        *pdf = p_h / (4 * dot(*wi, h));

        return MicrofacetBSDF::f(wo, *wi);
    }

//  *wi = cosineHemisphereSampler.get_sample(pdf);
//  return MicrofacetBSDF::f(wo, *wi);
}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Project 3-2: Part 1
  // Implement RefractionBSDF
    
    double eta_prime;
    if (wo.z >= 0) {
        eta_prime = 1.0/ior;
    } else {
        eta_prime = ior;
    }
    bool not_internal_refraction = refract(wo, wi, eta_prime);
    if (!not_internal_refraction) {
        return Vector3D();
    } else {
        *pdf = 1.0;
        return transmittance / abs_cos_theta(*wi) / (eta_prime * eta_prime);
    }
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Project 3-2: Part 1
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305
    double eta_prime;
    if (wo.z >= 0) {
        eta_prime = 1.0/ior;
    } else {
        eta_prime = ior;
    }

    bool not_internal_refraction = refract(wo, wi, eta_prime);
    if (!not_internal_refraction) {
        reflect(wo, wi);
        *pdf = 1.0;
        return reflectance / abs_cos_theta(*wi);
    } else {
        double cos_theta = abs(wo.z);
        double r0 = pow(((1 - ior) / (1 + ior)), 2);
        double r = r0 + (1-r0) * pow(1-cos_theta, 5);

        if (coin_flip(r)) {
            reflect(wo, wi);
            *pdf = r;
            return r * reflectance / abs_cos_theta(*wi);
        } else {
            refract(wo, wi, eta_prime);
            *pdf = 1-r;
            return (1-r) * transmittance / abs_cos_theta(*wi) / (eta_prime * eta_prime);
        }
    }

}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

  // TODO Project 3-2: Part 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
    Vector3D normal = Vector3D(0,0,1);
    *wi = -1 * wo + 2 * (wo * normal) * normal;


}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

  // TODO Project 3-2: Part 1
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
    double cos_theta = abs(wo.z);
    double cos_theta_prime = 1.0 - pow(ior, 2) * (1.0 - pow(cos_theta, 2));
    
    if (cos_theta_prime < 0) {
        return false;
    } else {
        double z_out;
        if (wo.z >= 0) {
            z_out = -1.0 * sqrt(1.0 - pow(ior, 2) * (1.0 - pow(wo.z, 2)));
        } else {
            z_out = 1.0 * sqrt(1.0 - pow(ior, 2) * (1.0 - pow(wo.z, 2)));
        }
        
        //store
        double x_out = -1.0 * ior * wo.x;
        double y_out = -1.0 * ior * wo.y;
        *wi = Vector3D(x_out, y_out, z_out);
        
        return true;
    }


}

} // namespace CGL
