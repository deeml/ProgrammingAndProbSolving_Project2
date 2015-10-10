package pb.g0;

import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;

import java.util.Random;
import java.util.ArrayList;
import java.util.Collections;

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

    private static final double G = Orbit.GM / Orbit.M;//6.67382967579392e-11;
    private static final double M = Orbit.M;//1.98855e30;

    // number of retries

    // print orbital information
    public void init(Asteroid[] asteroids, long time_limit)
    {
	if (Orbit.dt() != 24 * 60 * 60)
	    throw new IllegalStateException("Time quantum is not a day");
	this.time_limit = time_limit;
	prevNumAsteroids = asteroids.length;
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

	System.out.println("Angle: " + Double.toString(angle));
	System.out.println("Target angle: " + Double.toString(targetAngle));
	System.out.println("Speed: " + Double.toString(speed));
	System.out.println("Target speed: " + Double.toString(targetSpeed));

	/*System.out.println(vx);
	System.out.println(vy);
	System.out.println(target_vx);
	System.out.println(target_vy);
	System.out.println(push_vx);
	System.out.println(push_vy);
	System.out.println(push_angle);
	System.out.println(push_speed);*/
	//System.exit(0);
	ArrayList<Double> parameters = new ArrayList<Double>();
	parameters.add(new Double(push_speed));
	parameters.add(new Double(push_angle));
	return parameters;
    }

    // try to push asteroid
    public void play(Asteroid[] asteroids,
		     double[] energy, double[] direction)
    {
	time++;
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
	/*System.out.println("New iteration:");
	System.out.println(normv2);
	System.out.println( Math.sqrt( G * (M + outerAsteroid.mass) / r2));*/
	
	if (n < prevNumAsteroids) { //collision occurred, correct orbit of outermost asteroid
	    System.out.println("Collision!");
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
	    return;
	}
	
	if (collisionStarted) {return;}
	for (int i=n-2; i>=0; i--) {
	    int innerIndex = asteroids_ordered.get(i).index;	    
	    Asteroid innerAsteroid = asteroids[innerIndex];
	    if (innerAsteroid.orbit.a != innerAsteroid.orbit.b) {continue;}
	    double r1 = asteroids_ordered.get(i).getRadius();
	    Point v1 = innerAsteroid.orbit.velocityAt(time - innerAsteroid.epoch);
	    double normv1 = l2norm(v1);
	    double theta1 = Math.atan2(v1.y,v1.x);
	    double tH = Math.PI * Math.sqrt( Math.pow(r1+r2,3) / (8*Orbit.GM));
	    double thresh = innerAsteroid.radius() + outerAsteroid.radius();
	    // System.out.println(theta1 + " ---- " + theta2);
	    //System.out.println(  thresh / r2 - Math.abs(theta1 + Math.PI - theta2 - tH*omega2));
	    if ( Math.abs(theta1 + Math.PI - theta2 - tH*omega2) < thresh / r2) {
		double vnew = Math.sqrt(Orbit.GM / r1) * (Math.sqrt( 2*r2 / (r1+r2)) - 1);
		energy[innerIndex] = 0.5*asteroids[innerIndex].mass * vnew * vnew;
		direction[innerIndex] = theta1;
		collisionStarted = true;
	    }
	}
    }
}
