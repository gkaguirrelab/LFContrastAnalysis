cp = UnitCircleGenerate(30);

figure; hold on
scatter(cp(1,:),cp(2,:))

% parameters
p = [1,.5,45];

% Scale Matrix
S = diag(p(1:2));

% Rotation Matrix
V = deg2rotm(p(3))';

A_new = V' * S;
C =  A_new * cp;
C2 = inv(A_new) * cp;

[A_old,Ainv_old,Q] = EllipsoidMatricesGenerate(p,'dimension',2);
D =  Ainv_old * cp;
D2 = A_old * cp;

plot(C(1,:),C(2,:),'k')
plot(C2(1,:),C2(2,:),'k--')
plot(D(1,:),D(2,:),'r')
plot(D2(1,:),D2(2,:),'r--')
ylim([-2 2])
xlim([-2,2])
axis square

Q_new = A_new'*A_new;


r_old = diag(sqrt(cp'*Q*cp));
r_new = diag(sqrt(cp'*Q_new*cp));

figure; hold on
plot(r_old, 'r')
plot(r_new,'k')

