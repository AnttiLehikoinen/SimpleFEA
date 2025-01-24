classdef Mesh < handle
    properties
        nodes
        elements
    end
    methods
        function F = get_mapping_Jacobian_for_element(this, element_index)
            F = zeros(2,2);
            F(:,1) = this.nodes(:, this.elements(2, element_index)) - ...
                this.nodes(:, this.elements(1, element_index));
            F(:,2) = this.nodes(:, this.elements(3, element_index)) - ...
                this.nodes(:, this.elements(1, element_index));
        end

        function h = plot_elements(this, elements_to_plot, varargin)
            if isempty(elements_to_plot)
                elements_to_plot = 1:size(this.elements,2);
            end
            h = triplot(this.elements(:, elements_to_plot)', ...
                this.nodes(1,:), this.nodes(2,:), ...
                varargin{:});
        end

        function h = plot_nodes(this, nodes_to_plot, varargin)
            h = plot(this.nodes(1, nodes_to_plot), ...
                this.nodes(2, nodes_to_plot), varargin{:});
        end
    end
end