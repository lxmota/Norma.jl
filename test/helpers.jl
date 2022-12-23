function average_components(displacement::Vector{Float64})
    avg_x = mean(displacement[1:3:end])
    avg_y = mean(displacement[2:3:end])
    avg_z = mean(displacement[3:3:end])
    return [avg_x avg_y avg_z]    
end

function average_components(stress::Vector{Vector{Vector{Vector{Float64}}}})
    avg_xx = 0.0
    avg_yy = 0.0
    avg_zz = 0.0
    avg_yz = 0.0
    avg_xz = 0.0
    avg_xy = 0.0
    num_stress = 0
    for blk_index ∈ 1:length(stress)
        for blk_elem_index ∈ 1:length(stress[blk_index])
            for point ∈ 1:length(stress[blk_index][blk_elem_index])
                avg_xx += stress[blk_index][blk_elem_index][point][1]
                avg_yy += stress[blk_index][blk_elem_index][point][2]
                avg_zz += stress[blk_index][blk_elem_index][point][3]
                avg_yz += stress[blk_index][blk_elem_index][point][4]
                avg_xz += stress[blk_index][blk_elem_index][point][5]
                avg_xy += stress[blk_index][blk_elem_index][point][6]
                num_stress += 1
            end
        end
    end
    return [avg_xx avg_yy avg_zz avg_yz avg_xz avg_xy] ./ num_stress
end