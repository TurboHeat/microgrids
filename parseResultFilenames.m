function [parsed, files] = parseResultFilenames(resPath, resFormatSpec)
%{
EXAMPLE:
parsed = parseResultFolder();

figure(); histogram(parsed.Iteration, 'BinMethod', 'integers', 'EdgeColor', 'none').';

done = full(sparse(parsed.Iteration, parsed.("Fuel Price ID"), true, max(parsed.Iteration), max(parsed.("Fuel Price ID"))));
figure(); spy(done.'); ylabel('FuelId');

d2 = permute(reshape(done, [], 4, 6), [1,3,2]);
c = lines(4);
hF = figure('Position', [605,212,742,780]); hAx = gobjects(4,1); hT = tiledlayout(hF, 4, 1);
for id = 1:4
  hAx(id) = nexttile(hT);
  [x,y] = find( d2(:,:,id) );
  plot(hAx(id), x, y, '.', 'Color', c(id,:));
  title(hAx(id), "Building Type: " + string(BuildingType(id)));
end
xlim(hAx, [0 130]); xticks(hAx, 0:10:130);
ylim(hAx, [1 6]); yticks(hAx, 1:6);
xlabel(hAx, 'Case ID');
ylabel(hAx, 'Fuel Price ID');
grid(hAx, 'on');

%}
arguments
  resPath (1,1) string = "../Data/Results";
  resFormatSpec (1,1) string = "I%04u_AT%02d_AP%03u_B%1u_F%1u.mat";
end
% List files
files = dir(resPath);
files = struct2table(files(~[files.isdir]));

% Turn the list into a table
parsed = array2table(cell2mat(...
  cellfun(@(x)sscanf(x, resFormatSpec), files.name, 'UniformOutput', false).')...
  .', 'VariableNames', ["Iteration", "Algorithm Type", "Algorithm Parameters",...
  "Building Type", "Fuel Price ID"]);

end