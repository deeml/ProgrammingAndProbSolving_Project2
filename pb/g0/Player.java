package pb.g0;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Random;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Player implements pb.sim.Player {

    // used to pick asteroid and velocity boost randomly
    private Random random = new Random();

    // current time, time limit
    private long time = -1;
    private long time_limit = -1;

    // time until next push
    //private long time_of_collision = 0;
    private boolean collisionStarted = false;

    private int prevNumAsteroids = -1;

    private long accumulatorID;

    private static final double G = Orbit.GM / Orbit.M;//6.67382967579392e-11;
    private static final double M = Orbit.M;//1.98855e30;

    private int init_num_asteroids;
    private long time_gap; // Will look by this time into the future
    // number of retries

    HashMap<Long, push_move> push_queue = new HashMap<Long, push_move>();

    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit)
    {
	if (Orbit.dt() != 24 * 60 * 60)
	    throw new IllegalStateException("Time quantum is not a day");
	this.time_limit = time_limit;
	prevNumAsteroids = asteroids.length;
	init_num_asteroids = asteroids.length;
	int n = asteroids.length;
	ArrayList<asteroid_index> asteroids_ordered = new ArrayList <asteroid_index> ();
	Point astr_location = new Point();
	for (int i=0; i<n; i++) {
	    asteroids[i].orbit.positionAt(time - asteroids[i].epoch, astr_location);
	    asteroids_ordered.add(new asteroid_index(i, l2norm(astr_location)));
	}
	Collections.sort(asteroids_ordered);	
	accumulatorID = asteroids[asteroids_ordered.get(n-1).index].id;
    }

    private double distanceToSun(Asteroid x) {
	Point location = new Point();
	x.orbit.positionAt(time, location);
	return l2norm(location);
    }

    private double l2norm(Point p) {return Math.sqrt(p.x*p.x+p.y*p.y);}
    private double l2norm(double x, double y) {return Math.sqrt(x*x+y*y);}

    private ArrayList<Double> calculatePush(double speed, double angle, double targetSpeed, double targetAngle) {
	double vx = speed * Math.cos(angle);
	double vy = speed * Math.sin(angle);
	double target_vx = targetSpeed * Math.cos(targetAngle);
	double target_vy = targetSpeed * Math.sin(targetAngle);
	double push_vx = target_vx - vx;
	double push_vy = target_vy - vy;
	double push_angle = Math.atan2(push_vy, push_vx);
	double push_speed = l2norm(push_vx, push_vy);

	ArrayList<Double> parameters = new ArrayList<Double>();
	parameters.add(new Double(push_speed));
	parameters.add(new Double(push_angle));
	return parameters;
    }

    public boolean pathFree(Asteroid[] asteroids, ArrayList<asteroid_index> asteroids_ordered, Asteroid newAsteroid, int index, long startTime) {
	long t = startTime;
	int n = asteroids.length;
	Point p1 = new Point();
	Point p2 = new Point();
	while (t < time_limit) {
	newAsteroid.orbit.positionAt(t - newAsteroid.epoch, p1);
	    for (int i=index+1; i<n; i++) {
		asteroids[asteroids_ordered.get(i).index].orbit.positionAt(t - asteroids[asteroids_ordered.get(i).index].epoch, p2);
		if (l2norm(p1.x-p2.x,p1.y-p2.y) < asteroids[asteroids_ordered.get(i).index].radius() + newAsteroid.radius()) {
		    if (i < n-1) {
			return false;
		    }
		    else {
			return true;
		    }
		}
	    }
	    t++;
	}
	System.out.println("No collision...");
	return false;
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids,
		     double[] energy, double[] direction)
    {
	time++;
	if(push_queue.containsKey(time))
	{
		push_move p = push_queue.get(time);
		for (int i=0; i<asteroids.length; i++) {
		    if (p.id == asteroids[i].id) {
			energy[i] = p.energy;
			direction[i] = p.direction;
		    }
		}
	}

	int n = asteroids.length;
	ArrayList<asteroid_index> asteroids_ordered = new ArrayList();
	Point astr_location = new Point();
	for (int i=0; i<n; i++) {
	    asteroids[i].orbit.positionAt(time - asteroids[i].epoch, astr_location);
	    asteroids_ordered.add(new asteroid_index(i, l2norm(astr_location)));
	}
	Collections.sort(asteroids_ordered);

	if (n < prevNumAsteroids) { //collision occurred, correct orbit of outermost asteroid
	    boolean collided = true;
	    for (int i=0; i<asteroids.length; i++) {
		if (asteroids[i].id == accumulatorID) {
		    collided = false;
		    break;
		}
	    }
	    if (collided) {
		accumulatorID = -1;
		int accumulatorIndex = -1;
		for (int i=0; i<asteroids.length; i++) {
		    if (asteroids[i].id > accumulatorID) {
			accumulatorID = asteroids[i].id;
			accumulatorIndex = i;
		    }
		}		
		/*		Point location = new Point();
		Asteroid accumulator = asteroids[accumulatorIndex];
		accumulator.orbit.positionAt(time - accumulator.epoch, location);
		double r2 = l2norm(location); 
		double orbit_speed = Math.sqrt(Orbit.GM / r2);
		double tangent_theta = Math.PI/2 + Math.atan2(location.y, location.x);
		double omega2 = Math.sqrt(Orbit.GM / Math.pow(r2,3));
		Point v2 = accumulator.orbit.velocityAt(time - accumulator.epoch);
		double normv2 = l2norm(v2);
		double theta2 = Math.atan2(v2.y,v2.x);
		ArrayList<Double> parameters = calculatePush(normv2, theta2, orbit_speed, tangent_theta);
		energy[accumulatorIndex] = 0.5 * accumulator.mass * Math.pow(parameters.get(0),2);
		direction[accumulatorIndex] = parameters.get(1);*/
		prevNumAsteroids = n;
		collisionStarted = false;
		return;
	    }
	}
	
	if (collisionStarted) {return;}
	int accumulatorIndex = -1;
	for (int i=0; i<asteroids.length; i++) {
	    if (asteroids[i].id == accumulatorID) {
		accumulatorIndex = i;
	    }
	}

	if (accumulatorIndex == -1) {
	    System.out.println("No valid accumulator asteroid");
	    System.exit(0);
	}
	Asteroid accumulator = asteroids[accumulatorIndex];
	//	double omega2 = Math.sqrt(Orbit.GM / Math.pow(r2,3));
	
	double min_E = Double.MAX_VALUE;
	double min_dir = Double.MAX_VALUE;
	int min_index = -1;
	long min_time = 0;
	Point v1 = new Point();
	Point v2 = new Point();
	Point pmin = new Point();
	Point pmax = new Point();
	Point p = new Point();
	
	for (int i=n-1; i>=0; i--) {
	    int pushIndex = asteroids_ordered.get(i).index;	    
	    if (asteroids[pushIndex].id == accumulatorID) {continue;}
	    Asteroid pushAsteroid = asteroids[pushIndex];
	    double r1 = asteroids_ordered.get(i).getRadius();
	    time_gap = (2 * time_limit)/init_num_asteroids; //365*40
	    
	    for (long ft = 0 ; ft < time_gap ; ++ft) {
		long t = time + ft;
		double thresh = pushAsteroid.radius() + accumulator.radius();		
		accumulator.orbit.velocityAt(t - accumulator.epoch,v2);
		double theta2 = Math.atan2(v2.y,v2.x);
		double normv2 = l2norm(v2);
		theta2 = Math.atan2(v2.y,v2.x);		
		pushAsteroid.orbit.velocityAt(t - pushAsteroid.epoch,v1);
		double normv1 = l2norm(v1);
		double theta1 = Math.atan2(v1.y,v1.x);

		double rmin = Math.min(accumulator.orbit.a,accumulator.orbit.b);
		double rmax = Math.max(accumulator.orbit.a,accumulator.orbit.b);
	        long tHmin = (long) Math.round(Math.PI * Math.sqrt( Math.pow(r1+rmin,3) / (8*Orbit.GM)));
		long tHmax = (long) Math.round(Math.PI * Math.sqrt( Math.pow(r1+rmax,3) / (8*Orbit.GM)));

		accumulator.orbit.velocityAt(t+tHmin, pmin);
		accumulator.orbit.velocityAt(t+tHmax, pmax);

		double theta_min = (Math.atan2(pmin.y, pmin.x) + 2*Math.PI) % (2*Math.PI);
		double theta_max = (Math.atan2(pmax.y, pmax.x) + 2*Math.PI) % (2*Math.PI);	       

		if (accumulator.orbit.a == accumulator.orbit.b) {
		    double omega2 = Math.sqrt(Orbit.GM / Math.pow(rmin,3));
		    if (Math.abs(theta1+Math.PI - theta2-tHmin*omega2) < thresh / rmin) {
			double vnew = Math.abs(Math.sqrt(Orbit.GM / r1) * (Math.sqrt( 2*rmin / (r1+rmin)) - 1));
			double curr_E = 0.5*asteroids[pushIndex].mass*vnew*vnew;
			if  (curr_E < min_E){
			    min_E = curr_E;
			    min_index = pushIndex;
			    min_time = t;
			    if (r1 < rmin) { min_dir = theta1; }
			    else { min_dir = theta1+Math.PI;}
			}
			break;
		    }

		}
		
		else if ( (theta_min <= theta1 + Math.PI) && (theta1 + Math.PI <= theta_max) && ((theta_max - theta_min + 2*Math.PI) % (2*Math.PI)) < Math.PI) { // collision possible by varying r2 of orbit transfer
		    System.out.println("Collision possible with asteroid " + Integer.toString(i));
		    /*		    System.out.println(theta_min);
		    System.out.println(theta1+Math.PI);
		    System.out.println(theta_max);
		    System.out.println(rmin);
		    System.out.println(rmax);
		    System.out.println(pmin.y + "," + pmin.x); // 
		    System.out.println(pmax.y + "," + pmax.x);*/
		    double r2 = 0.5*(accumulator.orbit.a+accumulator.orbit.b);
		    long tH = (long)Math.round(Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM)));
		    accumulator.orbit.velocityAt(t+tH, p);
		    double theta_collide = (Math.atan2(p.y,p.x) + 2*Math.PI) % (2*Math.PI);
		    double step = 0.25 * Math.abs(accumulator.orbit.a - accumulator.orbit.b);		    
		    while ( Math.abs(theta_collide - theta1 - Math.PI) >= thresh / r2) { // search for colliding theta as a function of r2
			/*			System.out.println("New iteration");
			System.out.println(r2);
			System.out.println(theta_collide);
			System.out.println(theta1+Math.PI);*/
			if ( theta_collide < theta1 + Math.PI) {
			    r2 += step;
			}
			else {			    
			    r2 -= step;
			}
			tH = (long) Math.round(Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM)));
			accumulator.orbit.velocityAt(t+tH,p);
			theta_collide = (Math.atan2(p.y,p.x) + 2*Math.PI) % (2*Math.PI);
			step *= 0.5;
			//if (step < Math.sqrt(thresh)) {System.exit(0);}
		    }
		    double vnew = Math.abs(Math.sqrt(Orbit.GM / r1) * (Math.sqrt( 2*r2 / (r1+r2)) - 1));
		    double curr_E = 0.5*asteroids[pushIndex].mass*vnew*vnew;
		    if  (curr_E < min_E){
			min_E = curr_E;
			min_index = pushIndex;
			min_time = t;
			if (r1 < r2) { min_dir = theta1; }
			else { min_dir = theta1+Math.PI;}
		    }
		    break;
		}
	    }
	    
	}	      
	if(min_index != -1)
		{
		Asteroid a1 = null;
		try {
			a1 = Asteroid.push(asteroids[min_index], min_time, min_E, min_dir);
		}
		catch (pb.sim.InvalidOrbitException e) {
			System.out.println(" Invalid orbit: " + e.getMessage());
		}
		//if (pathFree(asteroids, asteroids_ordered, a1, min_index, min_time)){
		collisionStarted = true;
		push_queue.put(min_time, new push_move(min_index, asteroids[min_index].id, min_E, min_dir));
			//}
		System.out.println("Inserted into hashmap : "+push_queue.get(min_time).id + " | "+push_queue.get(min_time).energy + "|"+ push_queue.get(min_time).direction);
		}
    }
}
