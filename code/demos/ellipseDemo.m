% This is a demo for David to show why I am confused with how we are
% implementing the ellipse fit. 
%
% This will create a series of plots that made with the code currently used 
% in the QCM and plots with how I think the math should be. 

% first create a circle 
radius = 2;
x = radius*UnitCircleGenerate(30);

% Set the paramters of the ellipse. 
majorAxis = 1;
minorAxis = 0.5;
angle = 45;
p = [majorAxis,minorAxis,angle];

%% Current code - the function we use to create A,A_inv, and Q
[A_old,Ainv_old,Q] = EllipsoidMatricesGenerate(p,'dimension',2);
%{ 
This does:
    Scale matrix
    S = diag(p(1:2));

    Rotation Matrix - In the Brainard lab toolbox this has a transpose and I am not sure why
    since in the next step is is transposed again. 
    V = deg2rotm(p(3))';
    
    Create scale rotation matrix
    A = S*V'; I think this should be V'*S

    Create inverse scale rotation matrix
    Ainv = inv(A);

    Create the Q matrix 
    Q = A'*A;
%}
% transform the circle to an ellipse 
ep_old =  Ainv_old * x; % why doe we use the inverse?


%% The way that makes sense to me
%
% Create a scale matrix
S = diag(p(1:2));

% Create a rotation matrix
V = deg2rotm(p(3)); % no transpose

% Create scale rotation matrix
A_new = V * S;

% transform the circle to an ellipse 
ep_new =  A_new * x;



%% plot the results.
figure; 
hold on 
% plot the circle points 
scatter(x(1,:),x(2,:))

% plot the current way we do this 
plot(ep_old(1,:),ep_old(2,:),'r')

% plot the way that makes sense to me
plot(ep_new(1,:),ep_new(2,:),'k')

ylim([-4 4])
xlim([-4,4])
legend('circle points', 'current method', 'method that sense to michael','location','Northwest')
axis square

%% Look at ellipse distances 
% Current way
r_old = diag(sqrt(ep_old'*Q*ep_old));

% Way that makes sense to me
Q_new = A_new'*A_new;
r_new = diag(sqrt(ep_new'*Q_new*ep_new));

% plot it
figure; hold on
plot(r_old, 'r')
plot(r_new,'k--')

