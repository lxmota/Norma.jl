using LinearAlgebra
using Statistics
using Test
using Exodus

include("../src/Norma.jl")
include("helpers.jl")

@testset "transfer-operators" begin
    cp("../examples/contact/transfer_tests/cube_fine.g", "src_mesh.g", force=true)
    cp("../examples/contact/transfer_tests/cube_coarse.g", "dst_mesh.g", force=true)
    src_mesh = Exodus.ExodusDatabase("src_mesh.g", "r")
    dst_mesh = Exodus.ExodusDatabase("dst_mesh.g", "r")
    rm("src_mesh.g")
    rm("dst_mesh.g")
    #cube_fine contact side set 6, cube_coarse contact side set 5
    src_side_set_id = 5
    dst_side_set_id = 6
    src_T = get_boundary_traction_force(src_mesh, src_side_set_id)
    println("Number of nodes of the source side set ", length(src_T)) 
    dst_T_real = get_boundary_traction_force(dst_mesh, dst_side_set_id)
    println("Number of nodes of the destination side set ", length(dst_T_real)) 
    H = Norma.get_square_projection_matrix(src_mesh, src_side_set_id)
    L = Norma.get_rectangular_projection_matrix(dst_mesh, dst_side_set_id, src_mesh, src_side_set_id)
    dst_T = L * inv(H) * src_T
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
end