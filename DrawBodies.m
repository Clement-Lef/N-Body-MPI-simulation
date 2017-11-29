output = load('output.dat');

number_bodies = 40;


figure
for i = 1:2:2*number_bodies
    
    plot(output(:,i),output(:,i+1))
    hold on;
end
xlabel('Position x [arb.units]');
ylabel('Position y [arb.units]');
xlim([-2 2]);
ylim([-2 2]);