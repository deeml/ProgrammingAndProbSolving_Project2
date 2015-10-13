package pb.g0;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Random;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Player_checks implements pb.sim.Player {

    // used to pick asteroid and velocity boost randomly
    private Random random = new Random();

    // current time, time limit
    private long time = -1;
    private long time_limit = -1;

    // time until next push
    //private long time_of_collision = 0;
    private boolean collisionStarted = false;

    private int prevNumAsteroids = -1;

    private static final double G = Orbit.GM / Orbit.M;//6.67382967579392e-11;
    private static final double M = Orbit.M;//1.98855e30;

    private int init_num_asteroids;
    private long time_gap; // Will look by this time into the future
    // number of retries

    class push_move{
    	int index;
    	double energy;
    	double direction;
    	push_move(int index, double energy, double direction)
    	{
    		this.index = index;
    		this.energy = energy;
    		this.direction = direction;
    	}
    }

    HashMap<Long, push_move> push_queue = new HashMap<Long, push_move>();

    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit)
    {
	if (Orbit.dt() != 24 * 60 * 60)
	    throw new IllegalStateException("Time quantum is not a day");
	this.time_limit = time_limit;
	prevNumAsteroids = asteroids.length;
	init_num_asteroids = asteroids.length;
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
	System.out.println("  Year: " + (1 + time / 365));
	System.out.println("  Day: "  + (1 + time % 365));
	if(push_queue.containsKey(time))
	{

		// System.out.println("  ---- " );
		// System.out.println("  Year: " + (1 + time / 365));
		// System.out.println("  Day: "  + (1 + time % 365));
		// System.out.println("  ---- " );
		push_move p = push_queue.get(time);
		energy[p.index] = p.energy;
		direction[p.index] = p.direction; 
	}

	int n = asteroids.length;
	ArrayList<asteroid_index> asteroids_ordered = new ArrayList();
	Point astr_location = new Point();
	for (int i=0; i<n; i++) {
	    asteroids[i].orbit.positionAt(time - asteroids[i].epoch, astr_location);
	    asteroids_ordered.add(new asteroid_index(i, l2norm(astr_location)));
	}
	Collections.sort(asteroids_ordered);
	int outerIndex = asteroids_ordered.get(n-1).index;
	Asteroid outerAsteroid = asteroids[outerIndex];	
	double r2 = asteroids_ordered.get(n-1).getRadius();
	double omega2 = Math.sqrt(Orbit.GM / Math.pow(r2,3));
	Point v2 = outerAsteroid.orbit.velocityAt(time - outerAsteroid.epoch);
	double theta2 = Math.atan2(v2.y,v2.x);
	double normv2 = l2norm(v2);

	if (n < prevNumAsteroids) { //collision occurred, correct orbit of outermost asteroid
	    System.out.println("Collision at : "+"  Day: "  + (1 + time % 365)+"  Year: " + (1 + time / 365));
	    Point location = new Point();
	    outerAsteroid.orbit.positionAt(time - outerAsteroid.epoch, location);
	    double outerRadius = l2norm(location); // non circular orbit
	    //double orbit_speed = Math.sqrt( G*(M + outerAsteroid.mass) / outerRadius);
	    double orbit_speed = Math.sqrt(Orbit.GM / outerRadius);
	    //double orbit_speed = Math.sqrt( Orbit.GM * outerAsteroid.mass / outerRadius);
	    double tangent_theta = Math.PI/2 + Math.atan2(location.y, location.x);
	    ArrayList<Double> parameters = calculatePush(normv2, theta2, orbit_speed, tangent_theta);
	    energy[outerIndex] = 0.5 * outerAsteroid.mass * Math.pow(parameters.get(0),2);
	    direction[outerIndex] = parameters.get(1);
	    prevNumAsteroids = n;
	    collisionStarted = false;
	    // System.out.println("Energy for Orbit Correction: " + energy[outerIndex]);
	    // System.out.println("Mass being pushed: " + outerAsteroid.mass);
	    return;
	}
	
	if (collisionStarted) {return;}
	double min_E = Double.MAX_VALUE;
	double min_dir = Double.MAX_VALUE;
	int min_index = -1;
	long min_time = 0;
	// int count = 0;
	for (int i=n-2; i>=0; i--) {
	    int innerIndex = asteroids_ordered.get(i).index;	    
	    Asteroid innerAsteroid = asteroids[innerIndex];
	    if (innerAsteroid.orbit.a != innerAsteroid.orbit.b) {continue;} //if its not a circle, continue
	   	double r1 = asteroids_ordered.get(i).getRadius();
	   	time_gap = (2 * time_limit)/init_num_asteroids; //365*40
		for (long ft = 0 ; ft < time_gap ; ++ft) {
			long t = time + ft;
			Point v1 = innerAsteroid.orbit.velocityAt(t - innerAsteroid.epoch);
			double normv1 = l2norm(v1);
			double theta1 = Math.atan2(v1.y,v1.x);
			double tH = Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM));
			double thresh = innerAsteroid.radius() + outerAsteroid.radius();
			v2 = outerAsteroid.orbit.velocityAt(t - outerAsteroid.epoch);
			theta2 = Math.atan2(v2.y,v2.x);

			if ( Math.abs(theta1 + Math.PI - theta2 - tH*omega2) < thresh / r2) {
				// count ++;
				double vnew = Math.sqrt(Orbit.GM / r1) * (Math.sqrt( 2*r2 / (r1+r2)) - 1);
				double curr_E = 0.5*asteroids[innerIndex].mass*vnew*vnew;
				
				if  (curr_E < min_E){
				    min_E = curr_E;
					min_dir = theta1;
					min_index = innerIndex;
					min_time = t;
					}
				}
			}

		}	      

	if(min_index != -1)
		{
		//energy[min_index] = min_E;
		//direction[min_index] = min_dir;
		Asteroid a1 = null;
		try {
			a1 = Asteroid.push(asteroids[min_index], min_time, min_E, min_dir);
		}
		catch (pb.sim.InvalidOrbitException e) {
			System.out.println(" Invalid orbit: " + e.getMessage());
		}
		System.out.println("Pathfree output : "+pathFree(asteroids, asteroids_ordered, a1, min_index, min_time));
		if (pathFree(asteroids, asteroids_ordered, a1, min_index, min_time)){
			System.out.println("Inserted into hashmap : "+push_queue.get(min_time).index + " | "+push_queue.get(min_time).energy + "|"+ push_queue.get(min_time).direction);
			collisionStarted = true;
			push_queue.put(min_time, new push_move(min_index, min_E, min_dir));
		}
		// System.out.println("Inserted into hashmap : "+push_queue.get(min_time).index + " | "+push_queue.get(min_time).energy + "|"+ push_queue.get(min_time).direction);
		// System.out.println("  Year: " + (1 + min_time / 365));
		// System.out.println("  Day: "  + (1 + min_time % 365));
		// System.out.println("Energy for Push into Elliptical Orbit: " + energy[min_index]);
		}
    }
}
