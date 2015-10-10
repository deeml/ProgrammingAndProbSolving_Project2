package pb.g0;

public class CalendarDate{
	public long day;
	public long year;

	public CalendarDate()
	{
		this.day =  Long.valueOf(0);
		this.year = Long.valueOf(0);
	}

	public CalendarDate(long day, long year)
	{
		this.day = day;
		this.year = year;
	}
}