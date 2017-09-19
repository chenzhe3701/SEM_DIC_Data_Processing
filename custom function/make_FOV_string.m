% FOV = make_FOV_string(ri, rf, ci, cf, nDigits, sequence)
%
% chenzhe, 2017-05-30.
% make the string for FOV. [ri,rf,ci,cf] = start and finish index of r,c
% sequence = 'rc','snake',or 'raster'
function FOV = make_FOV_string(ri, rf, ci, cf, nDigits, sequence)

nR = rf - ri + 1;
nC = cf - ci + 1;
FOV = cell(nR,nC);

% the easy way is to generate, then just flip left and right.
% But I want the expressoin for the number of each position.
for iR = ri:rf
    for iC = ci:cf        
        switch sequence
            case {'RC', 'rc'}
                str = ['r',num2str(iR,['%.',num2str(nDigits),'d']),'c',num2str(iC,['%.',num2str(nDigits),'d'])];
            case {'Raster','raster'}
                num = (iR-ri)*nC + iC;
                str = num2str(num, ['%.',num2str(nDigits),'d']);
            case {'Snake', 'snake'}
                % if its odd row
                if mod((iR-ri),2)==1
                    % iC-ci = cf - iC_temp;
                    iC_temp = cf + ci - iC;
                else
                    iC_temp = iC;
                end
                num = (iR-ri)*nC + iC_temp;
                str = num2str(num, ['%.',num2str(nDigits),'d']);
            otherwise
                error('not supported sequence');
        end
        FOV{iR-ri+1, iC-ci+1} = str;
    end
end

end