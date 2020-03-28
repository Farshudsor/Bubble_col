results = load('results.mat');

sz = size(results.results);
hrs = 500;
Tsamp = 10;

Ce = zeros(sx(2), hrs);
Ca = zeros(sx(2), hrs);
Cz = zeros(sx(2), hrs);

D = zeros(sx(2), hrs/Tsamp);

T11 = zeros(sx(2), hrs/Tsamp);
T12 = zeros(sx(2), hrs/Tsamp);
T21 = zeros(sx(2), hrs/Tsamp);
T22 = zeros(sx(2), hrs/Tsamp);
T31 = zeros(sx(2), hrs/Tsamp);
T33 = zeros(sx(2), hrs/Tsamp);

Ce(:,:) = results.results.Ce;
Ca(:,:) = results.results.Ca;
Cx(:,:) = results.results.Cx;

D(:,:) = results.results.D;

T11(:,:) = results.results.T11;
T21(:,:) = results.results.T21;
T31(:,:) = results.results.T31;
T12(:,:) = results.results.T12;
T22(:,:) = results.results.T22;
T32(:,:) = results.results.T23;














