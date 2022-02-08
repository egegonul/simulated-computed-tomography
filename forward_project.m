clear all


load("square.mat");
load("SheppLogan.mat");
load("lena.mat");


img=lena;
img=img./max(img,[],'all');
M=size(img,1);


%specify step size and number of beams
num_beams=300;
step_size=1;
%create beams and angles
beams=linspace(-M/sqrt(2),M/sqrt(2),num_beams);
thetas=1:step_size:180;
projections=[];

all_data=cell(length(thetas),1);
for ww=1:length(thetas);
    theta=thetas(ww);
    projection=[];
    all_dist=cell(num_beams,1);
    all_adress=cell(num_beams,1);
    data=cell(num_beams,3);
    k=0;
    for t=beams
        k=k+1;
        %define x and y vectors
        x=(-M/2:1:M/2)';
        y=(-M/2:1:M/2)';
    
        %find intersections points for x and y
        yx=(t-x*cosd(theta))/(sind(theta)+0.00001);
        inds=yx<M/2&yx>-M/2;
        yx(abs(yx)<0.001)=0;
        yx=inds.*yx;
        inds=find(yx);
        relevant_x=[x(inds) yx(inds)];
        
        xy=(t-y*sind(theta))/(cosd(theta)+0.00001);
        xy(abs(xy)<0.001)=0;
        inds=xy<M/2&xy>-M/2;
        xy=inds.*xy;
        inds=find(xy);
        relevant_y=[xy(inds) y(inds)];
    
        %reduce the intersection points to relevant points
        %if ~isempty(relevant_x)&~isempty(relevant_y)
            relevant_points=unique([relevant_x;relevant_y],'rows');
            relevant_points=sortrows(relevant_points);
      % else
            %relevant_points=[];
       %end
        %calculate midpoints and distances only if there are nonzero relevant
        %points
        if(size(relevant_points,1)>1)
            %calculate distance and midpoint
            dist=[];
            mid=[];
            for i=1:size(relevant_points,1)-1
                x_cur=relevant_points(i,1);
                y_cur=relevant_points(i,2);
                x_next=relevant_points(i+1,1);
                y_next=relevant_points(i+1,2);
                dist=[dist; sqrt((x_next-x_cur)^2+(y_next-y_cur)^2)];
                mid=[mid;(x_cur+x_next)/2,(y_cur+y_next)/2];
            end
            data{k,2}=dist;
            
            
            %detect the adresses
            adress=[M/2-floor(mid(:,2)), M/2+ceil(mid(:,1))];
            data{k,3}=adress;
    
            %get the attenuation values 
            att=[];
            for i=1:size(adress,1)
                att=[att ;img(adress(i,1),adress(i,2))];
            end
        
            
            %compute the projection
            pj=att'*dist;
            data{k,1}=pj;
            projection=[projection; pj];
            
    
        else
            projection=[projection; 0];
        end
        
    end
    projections=[projections projection];
    all_data{ww}=data;
end


% data{1}=projection;
% data{2}=all_dist;
% data{3}=all_adress;
%plot the projection and the radon transform at the same angle together
subplot(2,1,1);
plot(beams,projection)
xlabel("t");
ylabel("p(t)");
title("projection at " +theta +" degrees using forward_project");
        
        
    subplot(2,1,2);
    plot(radon(img,theta))
    title("projection at " +theta +" degrees using Matlab's radon function")

save("projections.mat",'projections')


