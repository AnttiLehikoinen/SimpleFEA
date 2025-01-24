function [] = plot_flux_lines(a, msh, numberOfLines)
%plot_flux_density plots equipotential lines on a mesh.
% 
% plot_flux_density(a, msh, numberOfLines)
% plots a total of numberOfLines equipotential lines of the vector potential a 
% (defined on mesh msh) into the currently-selected figure.

N = numberOfLines;
Ne = size(msh.elements,2);

%potential values to be plotted
potentials = linspace(min(a), max(a), N+2);
potentials = potentials(2:(N+1));

hold on;
for kpot = 1:N
    pot = potentials(kpot);
    for elem = 1:Ne     
        %check if the equipotential lines goes through this element at all
        elementNodes = msh.elements(:, elem);
        if (max(a(elementNodes)) < pot) || (min(a(elementNodes)) > pot)
            %nope --> continue to the next element in list
            continue;
        end
        
        %drawing the equipotential line in the reference element:
            %denoting nodal values
            a1 = a(elementNodes(1)); a2 = a(elementNodes(2)); a3 = a(elementNodes(3));

            %finding where the equipotential line crosses the element boundaries,
            % ie. the the ksi-axis, eta-axis
            % and the 1-ksi line
            ksiis = [(pot-a1)/(a2-a1) 0 (pot-a3)/(a2-a3)];
            etas = [0 (pot-a1)/(a3-a1) 1-(pot-a3)/(a2-a3)];

            %finding the intersection points that are actually inside the
            %triangle
            I = find( (ksiis>=0).*(ksiis<=1).*(etas<=1).*(etas>=0) );
            ksiis = ksiis(I);
            etas = etas(I);

        %calculating mapping from reference element to global element
        %[B,b] = get_ElementwiseMapping(msh, elem);
        B = msh.get_mapping_Jacobian_for_element(elem);
        b = msh.nodes(:, msh.elements(1,elem));
        
        %calculating and plotting global equipotential lines
        globalCoordinates = B*[ksiis;etas] + [b b];
        plot( globalCoordinates(1,:)', globalCoordinates(2,:)', 'k')
    end
end

end

function [B,b] = get_ElementwiseMapping(msh, el)
%returns the mapping from the reference element to the global element el,
%such that g_global = B*x_ref + b

B = [msh.p(:,msh.t(2,el))-msh.p(:,msh.t(1,el)) msh.p(:,msh.t(3,el))-msh.p(:,msh.t(1,el))];
b = msh.p(:,msh.t(1,el));

end           
        