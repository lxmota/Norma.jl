#
#
function init(mesh::Mesh, domain::Domain)
    element_type = mesh.element_type
    num_int = mesh.num_int
    Na, dNa, w = isoparametric(element_type, num_int)
    domain.Na = Na
    domain.dNa = dNa
    domain.w = w
end
