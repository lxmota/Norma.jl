function [top_level_params, domain_data] = ...
    FullSchwarz(top_level_params, domain_data)

  num_domains = top_level_params.num_domains;
  num_steps = top_level_params.num_steps;
  step_interval = top_level_params.step_interval;
  num_stops = num_steps + 1;
  rel_tol_schwarz = top_level_params.rel_tol_schwarz;
  rel_tol_newton = top_level_params.rel_tol_newton;
  abs_tol_schwarz = top_level_params.abs_tol_schwarz;
  abs_tol_newton = top_level_params.abs_tol_newton;
  max_iter_schwarz = top_level_params.max_iter_schwarz;
  max_iter_newton = top_level_params.max_iter_newton;
  schwarz_interval = top_level_params.schwarz_interval;
  newton_interval = top_level_params.newton_interval;
  start_time = top_level_params.start_time;
  stop_time = top_level_params.stop_time;

  if top_level_params.problem_type == 2
    needs_mass = true;
  else
    needs_mass = false;
  end

  norms_displacements = zeros(num_domains, 1);
  norms_differences = zeros(num_domains, 1);
  norms_residual = zeros(num_domains, 1);
  time_diff = stop_time - start_time;
  dt = time_diff / num_steps;
  time_array = (start_time : dt : stop_time)';

  disp_errors = zeros(128, 1);

  %
  % Main loop
  %
  for stop = 1 : num_stops

    step = stop - 1;
    iteration_schwarz = 1;
    iterate_schwarz = true;
    display_step = mod(step, step_interval) == 0;

    if display_step == 1
      fprintf('Step: %d\n', step);
    end

    while iterate_schwarz == true

      display_schwarz = ...
          mod(iteration_schwarz, schwarz_interval) == 0 && display_step;

      % Go one domain at a time.
      for domain = 1 : num_domains

        this = domain_data{domain};
        num_elements = this.num_elements;
        ndof = this.ndof;
        num_dof = this.num_dof;
        num_int = this.num_int;
        dimension = this.dimension;
        nodes_per_element = this.nodes_per_element;

        free_index = this.free_index;

        displacements = this.displacements;
        prev_displacements = displacements;

        Na = this.Na;
        DNaDX = this.DNaDX;
        GradOp = this.GradOp;
        detDXDxi = this.detDXDxi;
        w = this.w;
        props = this.props;

        connectivity = this.connectivity;
        coordinates = this.coordinates;

        % set prescribed displacements here, before computing internal forces
        displacements = ...
            apply_bcs(domain_data, domain, stop, displacements, true);

        iteration_newton = 1;

        while true

          display_newton = mod(iteration_newton, newton_interval) == 0 ...
              && display_schwarz;

          % zero out the internal and external force vectors and stiffness
          % matrix
          Fint = zeros(num_dof, 1);
          Fext = zeros(num_dof, 1);
          K = zeros(num_dof, num_dof);

          if needs_mass == true
            M = zeros(num_dof, num_dof);
          end

          % initialize element stress/strain arrays
          e_stress = zeros(dimension, dimension, num_int);
          e_strain = zeros(dimension, dimension, num_int);
          Be = zeros(dimension, nodes_per_element, num_int);
          Beip = zeros(dimension, nodes_per_element);
          B = zeros(dimension * dimension, nodes_per_element * dimension);

          % loop through the elements and assemble their contributions to the
          % sttiffnes matrix and internal force vector          
          for i = 1 : num_elements

            % local gradient operator
            Be(:, :, :) = DNaDX(i, :, :, :);

            nodes = connectivity(i, :);
            X = coordinates(:, nodes)';
            
            % determine the degrees of freedom affected by this element 
            nd = ndof * nodes;
            nnm = [nd - 2; nd - 1; nd];
            dofs = reshape(nnm, [ndof * nodes_per_element, 1]);
            
            
            u = displacements(nnm');
            x = X + u;

            % loop over integration points
            nn = nodes_per_element * ndof;
            finte = zeros(nn, 1);
            ke = zeros(nn, nn);
            for l = 1 : num_int
              
              Beip(:, :) = Be(:, :, l);
              
              % compute deformation gradient and strain measure
              F = (Beip * x)';
              C = F' * F;
              
              % compute energy density, stress and moduli
              [~, S, CC] = Neohookean(props, C);
              
              P = F * S;
              AA = convect_tangent(CC, S, F);

              e_strain(:, :, l) = C;
              e_stress(:, :, l) = S;
              
              stress = reshape(P', dimension * dimension, 1);
              moduli = second_from_fourth(AA);
              B(:, :) = GradOp(i, l, :, :);
              
              Jw = detDXDxi(i, l) * w(l);

              % contribution to the internal force
              finte = finte + B' * stress * Jw;

              % contribution to the stiffness
              ke = ke + B' * moduli * B * Jw;

            end % integration points loop

            % place elem stress, state, strain into global array
            this.stress_history(:, :, i, :, stop) = e_stress;
            this.strain_history(:, :, i, :, stop) = e_strain;
            
            % assemble into global internal force vector
            Fint(dofs) = Fint(dofs) + finte;

            % if necessary, compute the element mass matrix
            if needs_mass == true
              rho = props.rho;
              me = zeros(nn, nn);
              if props.lumped == true
                meas = sum(detDXDxi(i, :));
                me = rho * meas * eye(nn, nn);
              else
                for k = 1 : num_int
                  for p = 1 : nodes_per_element
                    pp = ndof * (p - 1) + 1 : ndof * p;
                    for q = 1 : nodes_per_element
                      qq = ndof * (q - 1) + 1 : ndof * q;
                      mk = eye(ndof) * rho * w(k) ...
                           * detDXDxi(i, k) * Na(p, k) * Na(q, k);
                      me(pp, qq) = me(pp, qq) + mk;
                    end
                  end
                end
              end
            end

            % assemble into global stiffness
            K(dofs, dofs) = K(dofs, dofs) + ke;

            % assemble into global mass
            if needs_mass == true
              M(dofs, dofs) = M(dofs, dofs) + me;
            end

          end % elements loop

          % solve for the displacement increment
          F = Fint - Fext;
          R = F(free_index);
          delta = - K(free_index, free_index) \ R;
          displacements(free_index) = displacements(free_index) + delta;

          norm_delta = norm(delta);
          norm_disp = norm(displacements);

          if norm_disp > 0.0
            error = norm_delta / norm_disp;
          else
            error = norm_delta;
          end

          norm_residual = norm(R);
          norms_residual(domain) = norm_residual;
          norms_displacements(domain) = norm_disp;
          diff = displacements - prev_displacements;
          norms_differences(domain) = norm(diff);

          fmt_str = ['Subdomain: %d, Newton iter: %d, ' ...
                     '|\x03b4|=%0.17e, |\x03b4|/|u|=%0.17e, |r|=%0.17e\n'];
          if display_newton == 1
            fprintf(fmt_str, domain, iteration_newton, norm_delta, ...
                    error, norm_residual);
          end
          
          rel_pass = error <= rel_tol_newton;
          abs_pass = norm_delta <= abs_tol_newton;
          max_pass = iteration_newton >= max_iter_newton;

          if rel_pass || abs_pass || max_pass
            if display_schwarz == 1 && display_newton == 0
              fprintf(fmt_str, domain, iteration_newton, norm_delta, ...
                      error, norm_residual);
            end
            break;
          end

          iteration_newton = iteration_newton + 1;

        end % newton loop

        this.displacements = displacements;
        this.disp_history(:, stop) = displacements;

        domain_data{domain} = this;

      end % domains loop

      norm_displacements = norm(norms_displacements);
      norm_difference = norm(norms_differences);
      norm_residual = norm(norms_residual);

      if norm_displacements > 0.0
        error = norm_difference / norm_displacements;
      else
        error = norm_difference;
      end

      disp_errors(iteration_schwarz) = error;

      fmt_str = ['Step: %d, Schwarz iter: %d, |\x0394|=%0.17e, ' ...
      '|\x0394|/|U|=%0.17e, |R|=%0.17e'];

      if display_schwarz == 1
        fprintf(fmt_str, step, iteration_schwarz, norm_difference, ...
                error, norm_residual);
        for domain = 1 : num_domains
          fprintf(', |r_%d|=%0.17e', domain, norms_residual(domain));
        end
        fprintf('\n');
      end
      
      rel_pass = error <= rel_tol_schwarz;
      abs_pass = norm_difference <= abs_tol_schwarz;
      max_pass = iteration_schwarz >= max_iter_schwarz;
      
      if rel_pass || abs_pass || max_pass
        if display_step == 1 && display_schwarz == 0
          fprintf(fmt_str, step, iteration_schwarz, norm_difference, ...
                  error, norm_residual);
          for domain = 1 : num_domains
            fprintf(', |r_%d|=%0.17e', domain, norms_residual(domain));
          end
          fprintf('\n');
        end
        top_level_params.num_schwarz_iter = iteration_schwarz;
        break;
      end

      iteration_schwarz = iteration_schwarz + 1;

    end % iterate domains (Schwarz iteration)

  end % load step

  top_level_params.disp_errors = disp_errors;
  top_level_params.iteration_schwarz = iteration_schwarz;

  % write results
  basename = top_level_params.basename;
  for domain = 1 : num_domains
    this = domain_data{domain};
    mesh = this.mesh;
    num_elements = this.num_elements;
    mesh.Time = time_array;
    ndof = this.ndof;
    num_nodes = this.num_nodes;
    disp_history = this.disp_history;
    IZ = ndof * (1 : num_nodes)';
    IY = IZ - 1;
    IX = IZ - 2;
    disp_x = disp_history(IX, :);
    disp_y = disp_history(IY, :);
    disp_z = disp_history(IZ, :);
    mesh = AddNodalVar(mesh, 'disp_x', disp_x);
    mesh = AddNodalVar(mesh, 'disp_y', disp_y);
    mesh = AddNodalVar(mesh, 'disp_z', disp_z);
    
    stress = this.stress_history;
    
    stress_xx_1 = reshape(stress(1, 1, :, 1, :), num_elements, num_stops);
    stress_xy_1 = reshape(stress(1, 2, :, 1, :), num_elements, num_stops);
    stress_xz_1 = reshape(stress(1, 3, :, 1, :), num_elements, num_stops);
    
    stress_yx_1 = reshape(stress(2, 1, :, 1, :), num_elements, num_stops);
    stress_yy_1 = reshape(stress(2, 2, :, 1, :), num_elements, num_stops);
    stress_yz_1 = reshape(stress(2, 3, :, 1, :), num_elements, num_stops);

    stress_zx_1 = reshape(stress(3, 1, :, 1, :), num_elements, num_stops);
    stress_zy_1 = reshape(stress(3, 2, :, 1, :), num_elements, num_stops);
    stress_zz_1 = reshape(stress(3, 3, :, 1, :), num_elements, num_stops);

    mesh = AddElemVar(mesh, 'stress_xx_1', 1, stress_xx_1);
    mesh = AddElemVar(mesh, 'stress_xy_1', 1, stress_xy_1);
    mesh = AddElemVar(mesh, 'stress_xz_1', 1, stress_xz_1);

    mesh = AddElemVar(mesh, 'stress_yx_1', 1, stress_yx_1);
    mesh = AddElemVar(mesh, 'stress_yy_1', 1, stress_yy_1);
    mesh = AddElemVar(mesh, 'stress_yz_1', 1, stress_yz_1);
    
    mesh = AddElemVar(mesh, 'stress_zx_1', 1, stress_zx_1);
    mesh = AddElemVar(mesh, 'stress_zy_1', 1, stress_zy_1);
    mesh = AddElemVar(mesh, 'stress_zz_1', 1, stress_zz_1);
    
    domain_name = sprintf('%s_%02d', basename, domain - 1);
    output_name = strcat(domain_name, '.e');
    if exist(output_name, 'file') == 2
      delete(output_name);
    end
    exo_put(output_name, mesh);
  end

end