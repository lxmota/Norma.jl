using LinearAlgebra
using Statistics
using Test

include("../src/norma.jl")
include("helpers.jl")

@testset "transfer-traction" begin
    Exodus=Norma.exodus_module()
    src_mesh = Exodus.exodus("../examples/separate/transfer_tests/cube_fine.g")
    dst_mesh = Exodus.exodus("../examples/separate/transfer_tests/cube_coarse.g")
    #cube_fine contact side set 6, cube_coarse contact side set 5
    src_side_set_id = 6
    dst_side_set_id = 5
    src_T = get_boundary_traction_force(src_mesh, src_side_set_id)
    println("Number of nodes of the source side set ", length(src_T)) 
    dst_T_real = get_boundary_traction_force(dst_mesh, dst_side_set_id)
    println("Number of nodes of the destination side set ", length(dst_T_real)) 
    ratio = length(src_T) / length(dst_T_real)
    println("ratio ", ratio) 
    H = Norma.get_square_projection_matrix(src_mesh, src_side_set_id)
    dst_mesh.elem_to_blk_map = Norma.create_elem_to_blk_map(dst_mesh)
    L = Norma.get_rectangular_projection_matrix(src_mesh, src_side_set_id, dst_mesh, dst_side_set_id)
    dst_T = L * inv(H) * src_T
    #if ratio > 1.
    #    dst_mesh.elem_to_blk_map = Norma.create_elem_to_blk_map(dst_mesh)
    #    L = Norma.get_rectangular_projection_matrix(src_mesh, src_side_set_id, dst_mesh, dst_side_set_id)
    #    dst_T = L * inv(H) * src_T
    #else
    #    dst_mesh.elem_to_blk_map = Norma.create_elem_to_blk_map(src_mesh)
    #    L = Norma.get_rectangular_projection_matrix(dst_mesh, dst_side_set_id, src_mesh, src_side_set_id)
    #    dst_T =  transpose(L) * inv(H) * src_T
    #end    
    display([dst_T dst_T_real])
    rel_er = norm(dst_T - dst_T_real) / norm(dst_T_real)
    println("relative error ", rel_er) 
    @test norm(dst_T - dst_T_real) ≈ 0.0 atol = 1.0e-08
end