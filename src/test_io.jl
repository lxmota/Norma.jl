include("exodus.jl")
function test_io(input_file::String, output_file::String)
    rm(output_file, force = true)
    Exodus = exodus_module()
    mesh_struct = Exodus.copy_mesh(input_file, output_file)
    mesh_struct.put_time(1, 0.0)

    num_node_vars = mesh_struct.get_node_variable_number()
    disp_x_index = num_node_vars + 1
    disp_y_index = num_node_vars + 2
    disp_z_index = num_node_vars + 3
    num_node_vars += 3

    mesh_struct.set_node_variable_number(num_node_vars)

    mesh_struct.put_node_variable_name("disp_x", disp_x_index)
    mesh_struct.put_node_variable_name("disp_y", disp_y_index)
    mesh_struct.put_node_variable_name("disp_z", disp_z_index)
 
    num_element_vars = mesh_struct.get_element_variable_number()
    stress_xx_index = num_element_vars + 1
    stress_yy_index = num_element_vars + 2
    stress_zz_index = num_element_vars + 3
    stress_xy_index = num_element_vars + 4
    stress_xz_index = num_element_vars + 5
    stress_yz_index = num_element_vars + 6
    num_element_vars += 6

    mesh_struct.set_element_variable_number(num_element_vars)

    mesh_struct.put_element_variable_name("stress_xx", stress_xx_index)
    mesh_struct.put_element_variable_name("stress_yy", stress_yy_index)
    mesh_struct.put_element_variable_name("stress_zz", stress_zz_index)
    mesh_struct.put_element_variable_name("stress_xy", stress_xy_index)
    mesh_struct.put_element_variable_name("stress_xz", stress_xz_index)
    mesh_struct.put_element_variable_name("stress_yz", stress_yz_index)

    x, _, _ = mesh_struct.get_coords()
    num_nodes = length(x)
    disp_x = 1 .* ones(num_nodes)
    disp_y = 2 .* ones(num_nodes)
    disp_z = 3 .* ones(num_nodes)
    mesh_struct.put_variable_values("EX_NODAL", 1, "disp_x", 1, disp_x)
    mesh_struct.put_variable_values("EX_NODAL", 1, "disp_y", 1, disp_y)
    mesh_struct.put_variable_values("EX_NODAL", 1, "disp_z", 1, disp_z)

    elem_blk_ids = mesh_struct.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1 : num_blks
        blk_id = elem_blk_ids[blk_index]
        blk_conn = mesh_struct.get_elem_connectivity(blk_id)
        num_blk_elems = blk_conn[2]
        stress_xx = 1 .* ones(num_blk_elems)
        stress_yy = 2 .* ones(num_blk_elems)
        stress_zz = 3 .* ones(num_blk_elems)
        stress_xy = 4 .* ones(num_blk_elems)
        stress_xz = 5 .* ones(num_blk_elems)
        stress_yz = 6 .* ones(num_blk_elems)
        mesh_struct.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_xx", 1, stress_xx)
        mesh_struct.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_yy", 1, stress_yy)
        mesh_struct.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_zz", 1, stress_zz)
        mesh_struct.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_xy", 1, stress_xy)
        mesh_struct.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_xz", 1, stress_xz)
        mesh_struct.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_yz", 1, stress_yz)
    end

    mesh_struct.close()
end

input_file  = "/Users/amota/Repos/jlcm/examples/single/static-solid/unit-cube/cube.g"
output_file = "/Users/amota/Repos/jlcm/examples/single/static-solid/unit-cube/cube.e"

test_io(input_file, output_file)