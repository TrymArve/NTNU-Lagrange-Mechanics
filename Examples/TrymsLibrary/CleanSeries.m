%%% Clean Control Input Series

function[clean_U] = CleanSeries(U,tsim)

if length(U) > length(tsim)
i = 1;
U(:,1) = round(U(:,1),6);

while i < length(U)
ind = find(U(:,1) == U(i,1));
if length(ind) > 1
    U(ind(2:end),:) = [];
end
i = i + 1;
end

clean_U = interp1(U(:,1),U(:,2:end),tsim);
else
    clean_U = U;
end

end
