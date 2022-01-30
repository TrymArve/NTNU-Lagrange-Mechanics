%%% Clean Control Input Series



i = 1;
U(1,:) = round(U(1,:),6);

while i < length(U)
ind = find(U(1,:) == U(1,i));
if length(ind) > 1
    U(:,ind(2:end)) = [];
end
i = i + 1;
end

% disp('size U: '); size(U)
% disp('size tsim: '); size(tsim)
% disp('size U(1,:): '); size(U(1,:))
% disp('size U(2:end,:): '); size(U(2:end,:))

inputs = interp1(U(1,:)',U(2:end,:)',tsim);
% 
% disp('size inputs: '); size(inputs)

