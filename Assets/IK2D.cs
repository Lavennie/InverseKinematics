using System;
using System.Runtime.CompilerServices;
using UnityEngine;

public static class Helper
{
    /// <summary>
    /// 
    /// </summary>
    /// <param name="angle">In radians</param>
    /// <param name="radius"></param>
    /// <returns></returns>
    public static Vec2D ToVector(double angle, double radius)
    {
        return new Vec2D(Math.Cos(angle), Math.Sin(angle)) * radius;
    }
    public static double ToAngle(Vec2D point)
    {
        return Vec2D.SignedAngle(Vec2D.right, point);
    }
    public static IntervalD AngleInterval(Vec2D p1, Vec2D p2)
    {
        return new IntervalD(ToAngle(p1), ToAngle(p2));
    }
    public static Vec2D FollowChain(IKSegmentD[] segments, IndexRangeK fromToSplit, IntervalSide limitIndex1, IntervalSide limitIndex2)
    {
        Vec2D p1 = FollowChain(segments, new IndexRange(fromToSplit.startN, fromToSplit.k), limitIndex1, out double endAngle1);
        Vec2D p2 = FollowChain(segments, new IndexRange(fromToSplit.k, fromToSplit.endN), limitIndex2);
        return p1 + Vec2D.Rotate(p2, endAngle1);
    }
    public static Vec2D FollowChain(IKSegmentD[] segments, IndexRange fromTo, IntervalSide limitIndex)
    {
        return FollowChain(segments, fromTo, limitIndex, out double endAngle);
    }
    public static Vec2D FollowChain(IKSegmentD[] segments, IndexRange fromTo, IntervalSide limitIndex, out double endAngle)
    {
        endAngle = 0;
        Vector2 p = Vector2.zero;
        for (int i = fromTo.startN; i < fromTo.endN; i++)
        {
            if (limitIndex == IntervalSide.Min || limitIndex == IntervalSide.Max)
            {
                endAngle += segments[i].angleLimit[(int)limitIndex];
            }
            p += ToVector(endAngle, segments[i].length);
        }
        return p;
    }

    public static void DebugFollowChain(IKSegmentD[] segments, IndexRangeK fromToSplit, IntervalSide limitIndex1, IntervalSide limitIndex2, double startAngle, Color color)
    {
        Vec2D p1 = FollowChain(segments, new IndexRange(fromToSplit.startN, fromToSplit.k), limitIndex1, out double endAngle);
        DebugFollowChain(segments, new IndexRange(fromToSplit.startN, fromToSplit.k), limitIndex1, startAngle, Vec2D.zero, color, out double endAngle1);
        DebugFollowChain(segments, new IndexRange(fromToSplit.k, fromToSplit.endN), limitIndex2, endAngle1, p1, color);
    }
    public static void DebugFollowChain(IKSegmentD[] segments, IndexRange fromTo, IntervalSide limitIndex, double startAngle, Vec2D origin, Color color)
    {
        DebugFollowChain(segments, fromTo, limitIndex, startAngle, origin, color, out double endAngle);
    }
    public static void DebugFollowChain(IKSegmentD[] segments, IndexRange fromTo, IntervalSide limitIndex, double startAngle, Vec2D origin, Color color, out double endAngle)
    {
        endAngle = startAngle;
        Vec2D p = origin;
        for (int i = fromTo.startN; i < fromTo.endN; i++)
        {
            if (limitIndex == IntervalSide.Min || limitIndex == IntervalSide.Max)
            {
                endAngle += segments[i].angleLimit[(int)limitIndex];
            }
            Vec2D nextChain = ToVector(endAngle, segments[i].length);
            Debug.DrawLine(new Vector3((float)p.x, (float)p.y), new Vector3((float)(p.x + nextChain.x), (float)(p.y + nextChain.y)), color, 10);
            p += nextChain;
        }
    }

    public static CircleIntersectionD TwoCircleIntersection(CircleD c1, CircleD c2)
    {
        return TwoCircleIntersection(c1.center, c2.center, c1.radius, c2.radius);
    }
    public static CircleIntersectionD TwoCircleIntersection(CircleD c1, Vec2D c2, double r2)
    {
        return TwoCircleIntersection(c1.center, c2, c1.radius, r2);
    }
    public static CircleIntersectionD TwoCircleIntersection(Vec2D c1, Vec2D c2, double r1, double r2)
    {
        if (Vec2D.Distance(c1, c2) > r1 + r2 + 2 * Mathf.Epsilon)
        {
            return CircleIntersectionD.Miss;
        }
        if (c1.y != c2.y)
        {
            double k = -(c1.x - c2.x) / (c1.y - c2.y);
            double a = -((r1 * r1 - r2 * r2) - (c1.x * c1.x - c2.x * c2.x) - (c1.y * c1.y - c2.y * c2.y)) / (2 * (c1.y - c2.y));

            double d = ((c1.x + k * c1.y - k * a) * (c1.x + k * c1.y - k * a)) -
                (k * k + 1) * (c1.x * c1.x + c1.y * c1.y + a * a - 2 * a * c1.y - r1 * r1);
            // quadratic formula
            double x1 = ((c1.x + k * c1.y - k * a) + Math.Sqrt(d)) / (k * k + 1);
            double x2 = ((c1.x + k * c1.y - k * a) - Math.Sqrt(d)) / (k * k + 1);
            double y1 = k * x1 + a;
            double y2 = k * x2 + a;
            return new CircleIntersectionD(x1, y1, x2, y2);
        }
        else if (c1.x != c2.x)
        {
            // law of cosines without acos, since it gets canceled out when used in next line
            double d = Math.Abs(c1.x - c2.x);
            // TODO: swap epsilon for double variant
            if (Math.Abs(r1 + r2 - d) < (float)Mathf.Epsilon)
            {
                return new CircleIntersectionD((c1.x < c2.x) ? c1.x + r1 : c1.x - r1, c1.y, double.NaN, double.NaN);
            }
            else if (Math.Abs(r1 + d - r2) < Mathf.Epsilon)
            {
                return new CircleIntersectionD((c1.x < c2.x) ? c1.x - r1 : c1.x + r1, c1.y, double.NaN, double.NaN);
            }
            else if (Math.Abs(r2 + d - r1) < Mathf.Epsilon)
            {
                return new CircleIntersectionD((c1.x < c2.x) ? c1.x + r1 : c1.x - r1, c1.y, double.NaN, double.NaN);
            }
            double angle = (r1 * r1 + d * d - r2 * r2) / (2 * r1 * d);
            double offset = r1 * angle;
            double x = (c1.x < c2.x) ? c1.x + offset : c1.x - offset;
            CircleIntersectionD intersection = IntersectionLineCircle(90, new Vec2D(x - c1.x, 0), r1);
            intersection.p1 = new Vec2D(x, intersection.p1.y + c1.y);
            intersection.p2 = new Vec2D(x, intersection.p2.y + c1.y);
            return intersection;
            //return new CircleIntersection(new Vector2(x, c1.y), new Vector2(float.NaN, float.NaN));
        }
        else
        {
            return (r1 == r2) ? new CircleIntersectionD(c1.x, c1.y, r1, double.PositiveInfinity) : CircleIntersectionD.Miss;
        }
    }
    public static CircleIntersectionD IntersectionLineCircle(double angle, Vec2D center, double radius)
    {
        double radAngle = angle * Mathf.Deg2Rad;
        double cos = Math.Cos(radAngle);
        double sin = Math.Sin(radAngle);
        double y1 = sin * (center.y * sin + center.x * cos + Math.Sqrt(Math.Pow(center.y * sin + center.x * cos, 2) - center.x * center.x - center.y * center.y + radius * radius));
        double y2 = sin * (center.y * sin + center.x * cos - Math.Sqrt(Math.Pow(center.y * sin + center.x * cos, 2) - center.x * center.x - center.y * center.y + radius * radius));
        // TODO: sign depending on angle
        double x1 = center.x + Math.Sqrt(radius * radius - (y1 - center.y) * (y1 - center.y));
        double x2 = center.x - Math.Sqrt(radius * radius - (y2 - center.y) * (y2 - center.y));

        return new CircleIntersectionD(x1, y1, x2, y2);
    }
}


public class IK2D : MonoBehaviour
{

}


public struct CircleIntersectionD
{
    public enum Variant { Miss, Touching, Intersect, Covering }

    public Vec2D p1, p2;

    public CircleIntersectionD(double x1, double y1, double x2, double y2) : this(new Vec2D(x1, y1), new Vec2D(x2, y2)) { }
    public CircleIntersectionD(Vec2D p1, Vec2D p2)
    {
        this.p1 = p1;
        this.p2 = p2;
    }

    public Variant Type
    {
        get
        {
            if (double.IsNaN(p2.x))
            {
                if (double.IsNaN(p1.x))
                {
                    return Variant.Miss;
                }
                else
                {
                    return Variant.Touching;
                }
            }
            else if (double.IsInfinity(p2.y))
            {
                return Variant.Covering;
            }
            else
            {
                return Variant.Intersect;
            }
        }
    }
    public static CircleIntersectionD Miss { get { return new CircleIntersectionD(double.NaN, double.NaN, double.NaN, double.NaN); } }
}
public struct CircleD
{
    public Vec2D center;
    public double radius;

    public CircleD(Vec2D center, double radius)
    {
        this.center = center;
        this.radius = radius;
    }

    public override string ToString()
    {
        return $"({center}, {radius})";
    }
}
[System.Serializable]
public class IKSegmentD
{
    public double length;
    public IntervalD angleLimit;
    [Range(0, 1)] public double targetAngle = 0.5f;
    public bool animateTargetAngle = false;

    public IKSegmentD(double length, double angleMin, double angleMax)
    {
        this.length = length;
        this.angleLimit = new IntervalD(angleMin, angleMax);
    }
}

[System.Serializable]
public struct Vec2D
{
    public double x, y;

    public static readonly Vec2D zero = new Vec2D(0, 0);
    public static readonly Vec2D one = new Vec2D(1, 1);
    public static readonly Vec2D up = new Vec2D(0, 1);
    public static readonly Vec2D down = new Vec2D(0, -1);
    public static readonly Vec2D left = new Vec2D(-1, 0);
    public static readonly Vec2D right = new Vec2D(1, 0);


    public Vec2D(double x, double y)
    {
        this.x = x;
        this.y = y;
    }

    public static double Distance(Vec2D a, Vec2D b)
    {
        double num = a.x - b.x;
        double num2 = a.y - b.y;
        return Math.Sqrt(num * num + num2 * num2);
    }
    public static double Dot(Vec2D a, Vec2D b)
    {
        return a.x * b.x + a.y * b.y;
    }
    public static double Angle(Vec2D from, Vec2D to)
    {
        double num = Math.Sqrt(from.SqrMagnitude * to.SqrMagnitude);
        if (num < 1E-15f)
        {
            return 0f;
        }

        double num2 = Math.Clamp(Dot(from, to) / num, -1.0, 1.0);
        return Math.Acos(num2) * 57.29578;
    }
    public static double SignedAngle(Vec2D from, Vec2D to)
    {
        double num = Angle(from, to);
        double num2 = Math.Sign(from.x * to.y - from.y * to.x);
        return num * num2;
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="v"></param>
    /// <param name="angle">In radians</param>
    /// <returns></returns>
    public static Vec2D Rotate(Vec2D v, double angle)
    {
        double cos = Math.Cos(angle);
        double sin = Math.Sin(angle);
        return new Vec2D(v.x * cos - v.y * sin, v.x * sin + v.y * sin);
    }

    public double Magnitude { get { return Math.Sqrt(x * x + y * y); } }
    public double SqrMagnitude { get { return x * x + y * y; } }

    public static Vec2D operator +(Vec2D a, Vec2D b)
    {
        return new Vec2D(a.x + b.x, a.y + b.y);
    }
    public static Vec2D operator -(Vec2D a, Vec2D b)
    {
        return new Vec2D(a.x - b.x, a.y - b.y);
    }
    public static Vec2D operator *(Vec2D a, double x)
    {
        return new Vec2D(a.x * x, a.y * x);
    }
    public static Vec2D operator *(double x, Vec2D a)
    {
        return new Vec2D(a.x * x, a.y * x);
    }

}
[System.Serializable]
public struct IntervalD
{
    public double min, max;

    public IntervalD(double min, double max)
    {
        if (min < max)
        {
            this.min = min;
            this.max = max;
        }
        else
        {
            this.min = max;
            this.max = min;
        }
    }

    public static bool AreOverlapping(Interval i1, Interval i2)
    {
        // one end of i2 is inside i1
        return !Overlap(i1, i2).IsEmpty();
        //return i1.min == i2.min && i1.max == i2.max || i1.min < i2.min && i2.min < i1.max || i1.min < i2.max && i2.max < i1.max;
    }
    public static double ClampTo(double value, Interval interval)
    {
        return Math.Clamp(value, interval.min, interval.max);
    }
    public static Interval Overlap(Interval i1, Interval i2)
    {
        if (i1.min == i2.min && i1.max == i2.max)
        {
            return i1;
        }
        else if (i2.max < i1.min || i1.max < i2.min)
        {
            return Interval.Empty;
        }
        else
        {
            return new Interval(Math.Max(i1.min, i2.min), Math.Min(i1.max, i2.max));
        }
    }

    public bool IsEmpty()
    {
        return double.IsNaN(min) && double.IsNaN(max);
    }
    public bool IsPoint()
    {
        return min == max;
    }
    public bool Contains(double angle)
    {
        return min <= angle && angle <= max;
    }
    public double RandomInside()
    {
        System.Random rand = new System.Random();
        return min + rand.NextDouble() * (max - min);
    }

    public override string ToString()
    {
        return $"[{min}, {max}]";
    }

    public static Interval operator +(Interval i, double a)
    {
        return new Interval(i.min + a, i.max + a);
    }
    public static Interval operator +(double a, Interval i)
    {
        return new Interval(i.min + a, i.max + a);
    }

    public static Interval Empty { get { return new Interval(double.NaN, double.NaN); } }
    public double this[int i] { get { return (i % 2 == 0) ? min : max; } }
}