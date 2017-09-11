#
#
function bcDOF(mesh::Mesh, bc::BoundaryCondition)
    num_dof = mesh.dimension * mesh.num_nodes
    bc_dof = falses(num_dof, 1)
    num_nodesets = length(bc.nodesets)
    for nodeset_index = 1 : num_nodesets
        nodeset = bc.nodesets[nodeset_index]
        dof = bc.dof[nodeset_index]
        nodes = mesh.nodesets[nodeset]
        num_nodes = length(nodes)
        for node_index = 1 : num_nodes
            node = nodes[node_index]
            dof_index = (node - 1) * dimension + dof
            bc_dof[dof_index] = true
        end
    end
    return bc_dof
end

#
#
function indexDOF(mesh::Mesh)
    dbc_dof = bcDOF(mesh, mesh.dbc)
    sbc_dof = bcDOF(mesh, mesh.sbc)
    free_dof = .~(dbc_dof .! sbc_dof)
    return free_dof, dbc_dof, sbc_dof
end
