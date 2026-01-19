function S = assemble_stiffness_matrix(msh, reluctivity_per_element, shape_function_object)


N = size(msh.nodes, 2);
number_of_elements = size(msh.elements, 2);

quadrature_points = [0.33; 0.33];
quadrature_weights = 0.5;
number_of_quadrature_points = 1;

S = zeros(N, N);
for k_element = 1:number_of_elements
    for k_test = 1:3
        for k_shape = 1:3
            for k_quadrature = 1:number_of_quadrature_points
                %gradient of Ni as per lecture notes
                test_function_gradient_value = ...
                    shape_function_object.get_global_gradient_value(...
                    k_test, k_element, quadrature_points(:, k_quadrature));
                
                %gradient of Nj as per lecture notes
                shape_function_gradient_value = ...
                    shape_function_object.get_global_gradient_value(...
                    k_shape, k_element, quadrature_points(:, k_quadrature));
                
                %determinant of the Jacobian matrix of the mapping
                % from the reference into the global element k_element
                mapping_determinant = det(msh.get_mapping_Jacobian_for_element(k_element));

                %indices (row, column) into the stiffness matrix
                % as per lecture nodes, 'row' corresponds to i, and
                % 'column' corresponds to i
                row = msh.elements(k_test, k_element);
                column = msh.elements(k_shape, k_element);

                S(row, column) = S(row, column) + ...
                    quadrature_weights(k_quadrature) * ...
                    reluctivity_per_element(k_element) * ...
                    dot(test_function_gradient_value, shape_function_gradient_value) * ...
                    abs(mapping_determinant);
            end
        end
    end
end