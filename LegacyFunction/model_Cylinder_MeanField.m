function model = model_Cylinder_MeanField()
% return the mean-field model (Noack, 2011), which is a three state
% analytical model parameterized by Reynolds number (mu).

name = "Cylinder_MeanField";
nx = 3;

c = zeros(nx, 1);

mu = 1;
L = [mu, -1, 0;
    1, mu, 0;
    0, 0, -1];

Q = zeros(nx, nx, nx);
Q(3,1,1) = -0.5;
Q(1,3,1) = -0.5;
Q(2,3,2) = -0.5;
Q(3,2,2) = -0.5;
Q(:,:,3) = diag([1,1,0]);

model = class_Model_LosslessQuad(name, c, L, Q);
end