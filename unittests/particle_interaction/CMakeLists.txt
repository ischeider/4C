set(TESTNAME unittests_particle_interaction)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_particle_interaction_dem_contact_normal_test.cpp
    baci_particle_interaction_dem_contact_rolling_test.cpp
    baci_particle_interaction_dem_contact_tangential_test.cpp
    baci_particle_interaction_sph_artificialviscosity_test.cpp
    baci_particle_interaction_sph_density_correction_test.cpp
    baci_particle_interaction_sph_equationofstate_test.cpp
    baci_particle_interaction_sph_kernel_test.cpp
    baci_particle_interaction_sph_momentum_formulation_test.cpp
    baci_particle_interaction_utils_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} particle_interaction)
baci_target_requires_cyclic_legacy_linking(${TESTNAME})
