function [] = plot_flux_density(a, msh)
%plot_flux_density plots the flux density amplitude as color.
% 
% plot_flux_density(a, msh) plots the flux density amplitude, determined by
% the vector potential a and the mesh msh.

%getting the flux density in each element
Babs = getFluxDensity(a, msh);
Babs = repmat(Babs', 3, 1); %legacy fix for Matlab <2015

%describing the elements as polygons for Matlab
X = zeros(3, size(msh.elements,2));
Y = X;

for kn = 1:3
    X(kn,:) = msh.nodes(1, msh.elements(kn,:));
    Y(kn,:) = msh.nodes(2, msh.elements(kn,:));
end

%plotting
patch(X,Y, Babs, 'Linestyle', 'none');
colormap('jet');
colorbar;
axis tight;

end

function Babs = getFluxDensity(a, msh)
Ne = size(msh.elements, 2); %number of elements

Babs = zeros(Ne, 1);

shape_function = ShapeFunction();
shape_function.mesh = msh;
for ke = 1:Ne
    indices = msh.elements(:, ke);
    gradient_here = [0; 0];
    for k_shape = 1:3
        gradient_here = gradient_here + shape_function.get_global_gradient_value(...
            k_shape, ke, [0;0]) * ...
            a(indices(k_shape));
    end
    Babs(ke) = norm(gradient_here);    
end

end