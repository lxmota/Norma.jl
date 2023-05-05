using LinearAlgebra
using Statistics
using Test
using CPUTime

include("../src/norma.jl")
include("helpers.jl")

@testset "transfer-operators" begin
    Exodus=Norma.exodus_module()
    #src_mesh = Exodus.exodus("../examples/separate/transfer_tests/cube_coarse.g")
    #dst_mesh = Exodus.exodus("../examples/separate/transfer_tests/cube_fine.g")
    #src_mesh = Exodus.exodus("../examples/separate/transfer_tests/cube_fine_nonconf.g")
    #dst_mesh = Exodus.exodus("../examples/separate/transfer_tests/cube_coarse_nonconf.g")
    src_mesh = Exodus.exodus("../examples/separate/transfer_tests/cube_coarse_tet_1.g")
    dst_mesh = Exodus.exodus("../examples/separate/transfer_tests/cube_fine_tet_1-2.g")
    #cube_fine contact side set 6, cube_coarse contact side set 5
    src_side_set_id = 5
    dst_side_set_id = 6
    src_T = get_boundary_traction_force(src_mesh, src_side_set_id)
    #println(src_T)
    println("Number of nodes of the source side set ", length(src_T)) 
    dst_T_real = get_boundary_traction_force(dst_mesh, dst_side_set_id)
    println("Number of nodes of the destination side set ", length(dst_T_real)) 
    H = Norma.get_square_projection_matrix(src_mesh, src_side_set_id)
    src_mesh.elem_to_blk_map = Norma.create_elem_to_blk_map(src_mesh)
    L = Norma.get_rectangular_projection_matrix(dst_mesh, dst_side_set_id, src_mesh, src_side_set_id)
    dst_T = L * inv(H) * src_T
    #println(dst_T_real)
    display([dst_T dst_T_real])
    rel_er_tr = norm(dst_T - dst_T_real) / norm(dst_T_real)
    println("traction: relative error ", rel_er_tr) 
    @test norm(dst_T - dst_T_real) / norm(dst_T_real) ≈ 0.0 atol = 1.0e-08
    M = Norma.get_square_projection_matrix(dst_mesh, dst_side_set_id)
    src_u = ones(length(src_T))
    dst_u = inv(M) * L * src_u
    dst_u_real = ones(length(dst_T_real))
    display([dst_u dst_u_real])
    rel_er_disp = norm(dst_u - dst_u_real) / norm(dst_u_real)
    println("displacement: relative error ", rel_er_disp) 
    @test norm(dst_u - dst_u_real) / norm(dst_u_real) ≈ 0.0 atol = 1.0e-08
    #@time @CPUtime Norma.get_square_projection_matrix(src_mesh, src_side_set_id)
    #@time @CPUtime Norma.get_rectangular_projection_matrix(dst_mesh, dst_side_set_id, src_mesh, src_side_set_id)
end