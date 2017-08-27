function [top_level_params, domain_data] = Setup()

  top_level_params.basename = 'cuboid';

  %
  %    problem type: 1 - quasistatics
  %                  2 - dynamics
  %
  top_level_params.problem_type = 1;
  top_level_params.num_steps = 10;
  top_level_params.num_domains = 2;
  top_level_params.start_time = 0.0;
  top_level_params.stop_time = 1.0;

  top_level_params.rel_tol_schwarz = 1.0e-12;
  top_level_params.rel_tol_newton = 1.0e-12;
  top_level_params.abs_tol_schwarz =1.0e-12;
  top_level_params.abs_tol_newton = 1.0e-12;

  top_level_params.max_iter_schwarz = 1024;
  top_level_params.max_iter_newton = 1024;

  top_level_params.step_interval = 1;
  top_level_params.schwarz_interval = 1;
  top_level_params.newton_interval = 1;

  top_level_params.dimension = 3;
  top_level_params.ndof = 3;
  
  apply_initial_deformation = true;

  basename = top_level_params.basename;
  num_domains = top_level_params.num_domains;
  num_steps = top_level_params.num_steps;
  num_stops = num_steps + 1;

  domain_data = cell(num_domains);

  %
  % General parameters defined. Read each mesh and setup domain.
  %
  for domain = 1 : num_domains

    domain_name = sprintf('%s_%02d', basename, domain - 1);

    % Call each domain setup.
    fprintf('Reading mesh for subdomain %s ...\n', domain_name);
    fh = str2func(strcat('input_', domain_name));
    [mesh, props, bc] = fh(domain_name);

    fprintf('Initialization for subdomain %s ...\n', domain_name);

    num_nodes = size(mesh.Nodes.Coordinates, 1);
    dimension = size(mesh.Nodes.Coordinates, 2);
    num_elements = size(mesh.Blocks.Connectivity, 1);
    nodes_per_element = size(mesh.Blocks.Connectivity, 2);
    num_int = props.num_int;
    element_type = mesh.Blocks.ElementType;

    this.num_nodes = num_nodes;
    this.ndof = top_level_params.ndof;
    this.num_dof = num_nodes * top_level_params.ndof;
    this.dimension = dimension;
    this.element_type = element_type;
    this.num_elements = num_elements;
    this.nodes_per_element = nodes_per_element;
    this.num_int = num_int;
    this.mesh = mesh;
    this.props = props;
    this.bc = bc;
    
    num_dof = this.num_dof;

    [free_index, dbc_index, sbc_index, dbc_end_val] = input_dbcs(this);

    this.num_eqns = sum(free_index);

    dbc_val = zeros(size(dbc_end_val, 1), num_stops);

    for stop = 1 : num_stops
      dbc_val(dbc_index, stop) = ...
          (stop - 1) * dbc_end_val(dbc_index) / (num_stops - 1);
    end

    this.free_index = free_index;
    this.dbc_index = dbc_index;
    this.sbc_index = sbc_index;
    this.dbc_val = dbc_val;

    connectivity = mesh.Blocks.Connectivity;
    coordinates = mesh.Nodes.Coordinates';

    this.connectivity = connectivity;
    this.coordinates = coordinates;

    % Initialize arrays.
    this.displacements = zeros(num_dof, 1);
    this.disp_history = zeros(num_dof, num_stops);

    s2 = [dimension, dimension, num_elements, num_int, num_stops];
    
    target_size = props.target_size;
    this.target_size = target_size;
    
    this.Fini = InitialDeformation(num_stops, this, apply_initial_deformation);
    this.strain_history = zeros(s2);
    this.stress_history = zeros(s2);

    % initialize and compute integration rule weights, shape functions,
    % and local shape function derivatives
    [Na, DNaDxi, w] = isoparametric(element_type, num_int);

    this.Na = Na;
    this.DNaDxi = DNaDxi;
    this.w = w;

    % calculate gradient operator and jacobian determinant
    detDXDxi = zeros(num_elements, num_int);
    DNaDX = zeros(num_elements, dimension, nodes_per_element, num_int);
    a = dimension * dimension;
    b = nodes_per_element * dimension;
    GradOp = zeros(num_elements, num_int, a, b);

    for i = 1 : num_elements
      nodes = connectivity(i, :)';
      X = coordinates(:, nodes);
      for j = 1 : num_int
        DXDxiT = DNaDxi(:, :, j) * X';
        DNaDXij = DXDxiT \ DNaDxi(:, :, j);
        DNaDX(i, :, :, j) = DNaDXij;
        detDXDxi(i, j) = det(DXDxiT);
        GradOp(i, j, :, :) = gradient_operator(DNaDXij);
      end

    end

    this.GradOp = GradOp;
    this.DNaDX = DNaDX;
    this.detDXDxi = detDXDxi;

    domain_data{domain} = this;

  end % domain loop

  % Initialize Schwarz BCs
  for domain = 1 : num_domains
    domain_data = search_sbc(domain_data, domain);
  end

end % Setup function
