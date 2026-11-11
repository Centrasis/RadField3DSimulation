#pragma once
#include <string>
#include <vector>

namespace RadiationSimulation {
    namespace Analytical {

        /** Mass attenuation of the materials that are not water-like: for thin absorbers such as the image detector. */
        struct AttenuationTable {
            /** Whether a Geant4-derived table exists for the material (CsI, G4_WATER). */
            static bool has_table(const std::string& material);

            /** mu/rho in cm^2/g at e_kev (log-log interpolated), for materials with a table. */
            static float mass_attenuation(const std::string& material, float e_kev);

            /**
             * Optical thickness mu * thickness at perpendicular incidence, per field-spectrum bin (bin b spans
             * [b, b + 1) * bin_width_ev; the transmission is averaged over the bin). Materials without a table are
             * treated as water scaled by their density.
             */
            static std::vector<float> optical_thickness(const std::string& material, float density_g_cm3, float thickness_m,
                                                        int bins, float bin_width_ev);
        };
    }
}
