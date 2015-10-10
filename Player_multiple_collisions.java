package pb.g0;
import pb.sim.Point;
import pb.sim.Orbit;
import pb.sim.Asteroid;
import pb.sim.InvalidOrbitException;
// import java.util.Random;
import java.util.*;

public class Player implements pb.sim.Player {

	// used to pick asteroid and velocity boost randomly
	private Random random = new Random();

	// current time, time limit
	private long time = -1;
	private long time_limit = -1;

	// time until next push
	private long time_of_push = 0;

	// number of retries
	private int retries_per_turn = 1;
	private int turns_per_retry = 3;
	List<Integer> rem_asteroids = new ArrayList<>();
	CalendarDate[] collision_times = null;
	int new_asteroid = -1;
	int value = -1;
	Boolean collision = Boolean.FALSE;
	int count = -1;

	// print orbital information
	public void init(Asteroid[] asteroids, long time_limit)
	{
		collision_times = new CalendarDate[asteroids.length];
		Arrays.fill(collision_times, new CalendarDate());
		// Initially setting the indices of all the asteroids as 'rem_asteroids'//
		for (int i=0; i<asteroids.length; i++)
			rem_asteroids.add(i);
		Iterator<Integer> iter = rem_asteroids.iterator();
		System.out.println("Remaining asteroids are : --");
		while (iter.hasNext()) {
			System.out.println(iter.next());
		}
		System.out.println("--");
		if (Orbit.dt() != 24 * 60 * 60)
			throw new IllegalStateException("Time quantum is not a day");
		this.time_limit = time_limit;
	}

	// try to push asteroid
	public void play(Asteroid[] asteroids,
	                 double[] energy, double[] direction)
	{
		// if not yet time to push do nothing
		// if (++time <= time_of_push) return;

		++time;
		System.out.println("Year: " + (1 + time / 365));
		System.out.println("Day: "  + (1 + time % 365));
		CalendarDate current_date = new CalendarDate(1 + time % 365, 1 + time / 365);
		count = 0;
		for (int i=0; i<asteroids.length; i++)
		{	// add 1st guy not second. record second guy. decrement everyone above his value.
			System.out.println("Asteroid : "+ i + " | mass : "+ asteroids[i].mass);
			if ((current_date.day == collision_times[i].day) && (current_date.year == collision_times[i].year))
			{	//reset back asteroid collission time, if available asteroid value
				count ++;
				System.out.println("------------CD : " + current_date.day +"|"+current_date.year);
				collision_times[i] = new CalendarDate(); // resetting it

				if (count == 1)
				{// Collision flag has not been set yet; add the first guy and set it
					rem_asteroids.add(i); // Adding this asteroid back now
					System.out.println("Adding back : "+i);
					collision = Boolean.TRUE;
				}

				Iterator<Integer> iter = rem_asteroids.iterator();
				System.out.println("Remaining asteroids are : --");
				while (iter.hasNext()) {
					System.out.println(iter.next());
				}

				System.out.println("-- i --"+i);
				if(count == 2)
				{
					new_asteroid = i; // destroyed_asteroid //
					break;
				}
				
			}
		}
		if(collision){
			for (int a = 0; a<rem_asteroids.size(); a++)
			{
				System.out.println("a : " + a + " | " + rem_asteroids.get(a));
				// i was the value of the smaller asteroid
				if (rem_asteroids.get(a) > new_asteroid)
				{
					value = (int)rem_asteroids.get(a) - 1;
					rem_asteroids.set(a, value);
					System.out.println("new a : " + a + " | " + value);
				}
			}
			collision = Boolean.FALSE;
		}
		for (int retry = 1 ; retry <= retries_per_turn ; ++retry) {
			// pick a random asteroid and get its velocity
			// int i = random.nextInt(asteroids.length);
			System.out.println("REM AST : "+rem_asteroids.size());
			Iterator<Integer> iter1 = rem_asteroids.iterator();
			System.out.println("Remaining asteroids are : --");
			while (iter1.hasNext()) {
				System.out.println(iter1.next());
			}
			if(rem_asteroids.size() > 0)
			{
				int i = rem_asteroids.get(random.nextInt(rem_asteroids.size()));
				System.out.println("The asteroid picked is : "+i);
			
				for(int k=0; k<asteroids.length; k++)
					System.out.println("Asteroid : "+ k + " | mass : "+ asteroids[k].mass);
				Point v = asteroids[i].orbit.velocityAt(time);
				// add 5-50% of current velocity in magnitude
				System.out.println("Try: " + retry + " / " + retries_per_turn);
				double v1 = Math.sqrt(v.x * v.x + v.y * v.y);
				double v2 = v1 * (random.nextDouble() * 0.45 + 0.05);
				// System.out.println("  Speed: " + v1 + " +/- " + v2);
				// apply push at -π/8 to π/8 of current angle
				double d1 = Math.atan2(v.y, v.x);
				double d2 = d1 + (random.nextDouble() - 0.5) * Math.PI * 0.25;
				// System.out.println("  Angle: " + d1 + " -> " + d2);
				// compute energy
				double E = 0.5 * asteroids[i].mass * v2 * v2;
				// ---> try to push asteroid
				Asteroid a1 = null;
				try {
					a1 = Asteroid.push(asteroids[i], time, E, d2);
				} catch (InvalidOrbitException e) {
					System.out.println("  Invalid orbit: " + e.getMessage());
					continue;
				}
				// avoid allocating a new Point object for every position
				Point p1 = v, p2 = new Point();
				// search for collision with other asteroids
				// for (int j = 0 ; j != asteroids.length ; ++j) 
				Iterator<Integer> iterj = rem_asteroids.iterator();
				while (iterj.hasNext())
				{
					int j = iterj.next();
					if (i == j) continue;
					System.out.println("i : "+i+" | j : "+j);
					Asteroid a2 = asteroids[j];
					double r = a1.radius() + a2.radius();
					// look 10 years in the future for collision
					for (long ft = 0 ; ft != 3650 ; ++ft) {
						long t = time + ft;
						if (t >= time_limit) break;
						a1.orbit.positionAt(t - a1.epoch, p1);
						a2.orbit.positionAt(t - a2.epoch, p2);
						// if collision, return push to the simulator
						if (Point.distance(p1, p2) < r) {
							// Collision is going to happen between asteroid i & j!
							rem_asteroids.remove(new Integer(i));
							rem_asteroids.remove(new Integer(j));
							Iterator<Integer> iter = rem_asteroids.iterator();
							System.out.println("two asts are : " + i +"|" + j + "Remaining asteroids are : --");
							while (iter.hasNext()) {
								System.out.println(iter.next());
							}
							System.out.println("--");
							energy[i] = E;
							direction[i] = d2;
							// do not push again until collision happens
							time_of_push = t + 1;
							System.out.println("  Collision prediction between ! : " +i+ " & "+j );
							System.out.println("  Year: " + (1 + t / 365));
							System.out.println("  Day: "  + (1 + t % 365));
							// updating collision times array for asteroids i & j to value computed here;
							collision_times[i] = collision_times[j] = new CalendarDate((1 + t % 365), (1 + t / 365));
							System.out.println("----------coll value added : "+collision_times[i].day + " | "+collision_times[i].year);
							return;
						}
					}
				}
			}
			else
			{
				System.out.println("REM AST : "+rem_asteroids.size());
			}
			System.out.println("  No collision ...");
		}
		time_of_push = time + turns_per_retry;
	}
}
