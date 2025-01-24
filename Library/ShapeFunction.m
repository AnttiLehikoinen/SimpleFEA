classdef ShapeFunction < handle
    properties
        mesh
    end
    methods
        function dN = get_global_gradient_value(this, node_index, element_index, x_quad_local)
            J = this.mesh.get_mapping_Jacobian_for_element(element_index);

            dN_reference = this.get_reference_gradient_value(node_index, x_quad_local);
            dN = J' \ dN_reference;
        end

        function dN = get_reference_gradient_value(this, node_index, x_quad_local)
            if node_index == 1
                dN = [-1; -1];
            elseif node_index == 2
                dN = [1; 0];
            elseif node_index == 3
                dN = [0; 1];
            else
                error('Invalid node index.')
            end
        end

        function dN = get_global_value(this, node_index, element_index, x_quad_local)
            dN_reference = this.get_reference_value(node_index, x_quad_local);
            dN = dN_reference;
        end

        function N = get_reference_value(this, node_index, x_quad_local)
            if node_index == 1
                N = 1 - x_quad_local(1) - x_quad_local(2);
            elseif node_index == 2
                N = x_quad_local(1);
            elseif node_index == 3
                N = x_quad_local(2);
            else
                error('Invalid node index.')
            end
        end
    end
end