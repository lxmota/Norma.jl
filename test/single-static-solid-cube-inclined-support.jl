@testset "single-static-solid-cube-inclined-support" begin
    cp("../examples/single/static-solid/cube_compression_inclined_support_30_degrees/cube.yaml", "cube.yaml", force=true)
    cp("../examples/single/static-solid/cube_compression_inclined_support_30_degrees/cube.g", "cube.g", force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml")
    rm("cube.g")
    rm("cube.e")
    
    # These displacements are obtained from an identical problem with
    # BCs coincident with the lab basis
    reference_displacements = [
        0.19999999999999996, 0.0, 0.0, 0.19999999999999996, 
        0.0, 0.0, 0.19999999999999996, 
        0.0, 0.0, 0.19999999999999996, 0.0, 0.0, 
        -0.002490510706379005, -0.054808932034627424, 
        0.054808932034627285, 0.008794796508373306, 
        -0.05449689589044056, -2.8756551738792074e-17, 
        0.020548157939405773, -7.567237356520244e-17,
        -4.131059890383733e-18, 0.008794796508373347, 
        -4.2461283759579046e-17, 0.05449689589044062, 
        0.19999999999999996, 0.0, 0.0, 
        0.19999999999999996, 0.0, 0.0, 
        -0.0024905107063790998, -0.05480893203462743, 
        -0.05480893203462725, 0.008794796508373271,
        -4.832171592315198e-17, -0.054496895890440634, 
        0.19999999999999996, 0.0, 0.0, 0.19999999999999996, 
        0.0, 0.0, 0.008794796508373334, 
        0.05449689589044063, -5.145650465803633e-17,
        -0.002490510706379087, 0.05480893203462731, 
        0.0548089320346273, 0.19999999999999996, 0.0,
        0.0, -0.002490510706378977, 
        0.05480893203462727, -0.05480893203462724, -0.2,
        -0.052612760978146725, 
        0.052612760978146586, -0.2, -0.05432895977117359,
        9.641660437651523e-18, -0.2, 
        -3.8460049041556104e-17, -3.566854429999362e-17, 
        -0.2, -7.231031457420035e-17, 
        0.05432895977117351, -0.2, -0.05261276097814662,
        -0.052612760978146676, -0.2, 
        -1.0483851671759644e-16, -0.05432895977117339,
        -0.2, 0.054328959771173434, 
        -2.6230761947053614e-17, -0.2, 0.052612760978146704,
        0.05261276097814651, -0.2, 
        0.05261276097814652, -0.05261276097814667
    ]

    # Rotate these displacements
    angle = 30. * π / 180
    c = cos(angle)
    s = sin(angle)
    # Rotate about z
    local_rotation_matrix = [ c s 0; -s c 0 ; 0 0 1]
    global_rotation = zeros((81,81))
    for i in range(1,27)
        base = (i-1)*(3) + 1
        global_rotation[base:(base+2), base:(base+2)] = local_rotation_matrix
    end

    correct_displacements = global_rotation' * reference_displacements

    local_stiffness_displacements = integrator.displacement
    # Displacements in the integrator are stored in the local frame
    # (i.e., inclined support frame). Use the model global_transform
    # to return them to the global frame
    displacements = model.global_transform' * integrator.displacement

    # Assert the displacement array matches the reference displacements
    @test displacements ≈ correct_displacements atol=1e-6
    
end
