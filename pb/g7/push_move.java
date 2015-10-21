package pb.g7;

public class push_move{
    long id;
    int index;
    double energy;
    double direction;
    double density;
    long time;
    double radius;
    public push_move(int index, long id, double energy, double direction)
    {
	this.index = index;
	this.id = id;
	this.energy = energy;
	this.direction = direction;
    }
    public push_move(int index, long id, double energy, double direction, long time, double density, double radius)
    {
	this.index = index;
	this.id = id;
	this.energy = energy;
	this.direction = direction;
	this.time = time;
	this.density = density;
	this.radius = radius;
    }
}
