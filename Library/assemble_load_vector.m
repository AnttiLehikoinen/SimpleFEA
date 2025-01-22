function F = assemble_load_vector(msh, current_density, shape_function)


N = size(msh.nodes, 2);
number_of_elements = size(msh.elements, 2);

quadrature_points = [1/6 1/6;1/6 2/3;2/3 1/6]';
quadrature_weights = [1/6 1/6 1/6];
number_of_quadrature_points = 3;

F = zeros(N, 1);
for k_element = 1:number_of_elements
    for k_test = 1:3
        for k_quadrature = 1:number_of_quadrature_points
            test_function_value = ...
                shape_function.get_global_value(...
                k_test, k_element, quadrature_points(:, k_quadrature));
            mapping_determinant = det(msh.get_mapping_matrix_for_element(k_element));

            row = msh.elements(k_test, k_element);

            F(row) = F(row) + ...
                quadrature_weights(k_quadrature) * ...
                current_density(k_element) * ...
                test_function_value * ...
                abs(mapping_determinant);
        end
    end
end

end