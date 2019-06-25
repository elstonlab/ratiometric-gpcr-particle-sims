function analyze_screen()

activity_mode={'ratiometric','nonratiometric'};
densityGradient_modes = {'align','oppose','orthogonal','uniform','diagonal'};


for j=1:numel(activity_mode)
    for k=1:50
        for dg=1:5
            filebase=sprintf('20190328_R10000_%s_seed%i_denGrad%s',activity_mode{j},k,densityGradient_modes{dg});
            disp(filebase);
            [vec{1,j,k,dg},...
             t{1,j,k,dg},...
             nGa{1,j,k,dg},...
             confidence{1,j,k,dg},...
             norm_proj{1,j,k,dg},...
             xyzAlignment{1,j,k,dg},...
             xyAlignment{1,j,k,dg}]=quantify_v2(sprintf('/Volumes/SATAdrive/gpcr/20190330/N10000/%s.xyz',filebase),62,'Ga',2500); 
        end
    end
end

save('analyzed_data');
generate_csv_angleOverTime(t,xyAlignment,'angles_over_time')

end

function generate_csv_angleOverTime(t,xyAlignment,csvbasename)
% creates a csv file containing the angle-over-time data across all realizations,
% for each class of conditions.
dg_modes={'align','oppose','orthogonal','uniform','diagonal'};
modeltypes={'ratiometric','nonratiometric'};

for dg=1:5
    for j=1:2
        xyAlignmentmat = cell2mat(cellfun(@(x) x,squeeze(xyAlignment(1,j,:,dg))','uniformoutput',false));
        data_to_write = [t{1,j,1,dg}(2:end),xyAlignmentmat(2:end,:)];
        csvwrite(sprintf('%s_%s_%s.csv',csvbasename,dg_modes{dg},modeltypes{j}),data_to_write);
    end
end

end
