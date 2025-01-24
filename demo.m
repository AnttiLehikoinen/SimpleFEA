%demo Super-simple 2D FEA implementation
%
% Super-inefficient in all ways, too.
%
% Analyses (in 2D magnetostatics) a rectangular bus bar carrying some
% current.

%adding library functions
addpath('Library')

%loading mesh data
data = load('mesh_data.mat');

%initializing mesh
msh = Mesh();
msh.nodes = data.nodes;
msh.elements = data.elements;

%indices of elements and nodes
dirichlet_nodes = data.Dirichlet_nodes;
bar_elements = data.bar_elements;
air_elements = data.air_elements;

%sizes
number_of_elements = size(msh.elements, 2);
number_of_nodes = size(msh.nodes, 2);

%elementwise reluctivity and current density
mu0 = pi*4e-7; %vacuum permeability
reluctivity = ones(1, number_of_elements) / mu0;
current_density = zeros(1, number_of_elements);
current_density(bar_elements) = 6e6; %6 A per mmsq

%plotting
figure(1); clf; hold on; box on; axis equal;
hbar = msh.plot_elements(bar_elements, 'r');
hair = msh.plot_elements(air_elements, 'b');
hbnd = msh.plot_nodes(dirichlet_nodes, 'md', 'markersize', 10);
legend([hbar, hair, hbnd], 'Busbar', 'Air', 'Dirichlet boundary');


%assembling matrices
shape_function_object = ShapeFunction();
shape_function_object.mesh = msh;
S = assemble_stiffness_matrix(msh, reluctivity, shape_function_object);
S = sparse(S); %for faster solution of the linear problem
F = assemble_load_vector(msh, current_density, shape_function_object);

%solving free variables
free_nodes = setdiff(1:number_of_nodes, dirichlet_nodes);
A_free = S(free_nodes, free_nodes) \ F(free_nodes);

%parsing entire solution
A = zeros(number_of_nodes, 1);
A(free_nodes) = A_free;

%plotting solution
figure(2); clf; hold on; box on;
trimesh(msh.elements', msh.nodes(1,:), msh.nodes(2,:), A);

%plotting flux lines and flux density
figure(3); clf; hold on; box on; axis equal tight;
title('Flux lines and flux density')
plot_flux_density(A, msh);
plot_flux_lines(A, msh, 10);




