using UnityEngine;
using System.Collections.Generic;
using System.Collections;
#if UNITY_EDITOR
using UnityEditor;
#endif

public enum IntervalSide : byte { Min = 0, Max = 1 }

public static class IKUtility
{
    public static Vector2 ToVector(float angle, float radius)
    {
        return new Vector2(Mathf.Cos(angle * Mathf.Deg2Rad), Mathf.Sin(angle * Mathf.Deg2Rad)) * radius;
    }
    public static float ToAngle(Vector2 point)
    {
        return Vector2.SignedAngle(Vector2.right, point);
    }
    public static Interval AngleInterval(Vector2 p1, Vector2 p2)
    {
        return new Interval(ToAngle(p1), ToAngle(p2));
    }
    public static Vector2 FollowChain(IKSegment[] segments, int k, IntervalSide limitIndex)
    {
        return FollowChain(segments, k, k, limitIndex);
    }
    public static Vector2 FollowChain(IKSegment[] segments, int k, int n, IntervalSide limitIndex)
    {
        Vector2 p = Vector2.zero;
        float angleSum = 0;
        for (int i = 0; i < k; i++)
        {
            angleSum += segments[i].angleLimit[(int)limitIndex];
            p += ToVector(angleSum, segments[i].length);
        }

        if (k < n)
        {
            angleSum += segments[k].angleLimit[(int)limitIndex];
            for (int i = k; i < n; i++)
            {
                p += ToVector(angleSum, segments[i].length);
            }
        }
        return p;
    }

    public static Vector2 TwoCircleIntersection(Vector2 c1, Vector2 c2, float r1, float r2)
    {
        if (Vector2.Distance(c1, c2) > r1 + r2)
        {
            return new Vector2(float.NaN, float.NaN);
        }
        if (c1.y != c2.y)
        {
            float k = -(c1.x - c2.x) / (c1.y - c2.y);
            float a = -((r1 * r1 - r2 * r2) - (c1.x * c1.x - c2.x * c2.x) - (c1.y * c1.y - c2.y * c2.y)) / (2 * (c1.y - c2.y));

            float d = ((c1.x + k * c1.y - k * a) * (c1.x + k * c1.y - k * a)) -
                (k * k + 1) * (c1.x * c1.x + c1.y * c1.y + a * a - 2 * a * c1.y - r1 * r1);
            // quadratic formula
            float x1 = ((c1.x + k * c1.y - k * a) + Mathf.Sqrt(d)) / (k * k + 1);
            float x2 = ((c1.x + k * c1.y - k * a) - Mathf.Sqrt(d)) / (k * k + 1);
            float y1 = k * x1 + a;
            float y2 = k * x2 + a;
            return new Vector2(x1, y1);
        }
        // TODO: when line is vertical (c1.y == c2.y)
        return new Vector2(float.NaN, float.NaN);
    }
    public static Vector2 IntersectionLineCircle(float angle, Vector2 center, float radius)
    {
        float radAngle = angle * Mathf.Deg2Rad;
        float cos = Mathf.Cos(radAngle);
        float sin = Mathf.Sin(radAngle);
        float y;
        y = sin * (center.y * sin + center.x * cos + Mathf.Sqrt(Mathf.Pow(center.y * sin + center.x * cos, 2) - center.x * center.x - center.y * center.y + radius * radius));
        float x;
        x = center.x + Mathf.Sqrt(radius * radius - (y - center.y) * (y - center.y));

        return new Vector2(x, y);
    }
}

public class IKPolar : MonoBehaviour
{
    public bool debugChainBend = true;
    public bool debugTargetLine = true;

    public IKSegment[] segments;
    private Reach[] data;

    private Vector2 mousePos = Vector2.zero;
    private Vector2 reach = Vector2.zero;

    private void Awake()
    {
        Precalculate();
    }

    private void Update()
    {
        if (debugChainBend)
        {
            // draw min/max bend
            Vector2 pMin = Vector2.zero;
            Vector2 pMax = Vector2.zero;
            Interval angleSum = new Interval(0, 0);
            for (int i = 0; i < segments.Length; i++)
            {
                angleSum.min += segments[i].angleLimit[0];
                angleSum.max += segments[i].angleLimit[1];
                Vector2 dirMin = new Vector2(Mathf.Cos(angleSum.min * Mathf.Deg2Rad), Mathf.Sin(angleSum.min * Mathf.Deg2Rad)) * segments[i].length;
                Vector2 dirMax = new Vector2(Mathf.Cos(angleSum.max * Mathf.Deg2Rad), Mathf.Sin(angleSum.max * Mathf.Deg2Rad)) * segments[i].length;
                Debug.DrawRay(pMin, dirMin, Color.green, 0.01f);
                Debug.DrawRay(pMax, dirMax, Color.green, 0.01f);
                pMin += dirMin;
                pMax += dirMax;
            }
        }

        mousePos = Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, -Camera.main.transform.position.z));
        float angle = Mathf.Asin(mousePos.y / mousePos.magnitude) * Mathf.Rad2Deg;
        mousePos = ClampTarget(mousePos);


        Interval reachInterval = data[0].Get(segments, angle);
        reach = new Vector2(reachInterval.min, reachInterval.max);
        Debug.DrawRay(mousePos.normalized * reach.x, mousePos.normalized * Mathf.Max(reach.y - reach.x, Mathf.Epsilon), Color.red, 0.01f);
        if (debugTargetLine)
        {
            Debug.DrawRay(Vector3.zero, mousePos, Color.yellow, 0.01f);
        }

        //Align(mousePos);
    }
#if UNITY_EDITOR
    private void OnGUI()
    {
        EditorGUI.LabelField(new Rect(0, 0, 300, EditorGUIUtility.singleLineHeight), 
            $"Mouse position: {mousePos}");
        EditorGUI.LabelField(new Rect(0, EditorGUIUtility.singleLineHeight, 300, EditorGUIUtility.singleLineHeight), 
            $"Reach: {reach}");
    }
#endif

    public void Precalculate()
    {
        data = new Reach[Mathf.CeilToInt(segments.Length / 2.0f)];
        int i = 0;
        for (int k = segments.Length; k > 0; k -= 2)
        {
            data[i] = new Reach(segments, k);
            i++;
        }
    }

    public Vector2 ClampTarget(Vector2 target)
    {
        float angle = IKUtility.ToAngle(target);
        float radius = target.magnitude;

        angle = Interval.ClampTo(angle, data[0].ValidInterval);
        radius = Interval.ClampTo(radius, data[0].Get(segments, angle));
        return radius * new Vector2(Mathf.Cos(angle * Mathf.Deg2Rad), Mathf.Sin(angle * Mathf.Deg2Rad));
    }
    public bool Align(Vector3 target)
    {
        return false;
    }
}

[System.Serializable]
public class Reach
{
    private MinReachData minData;
    private MaxReachData maxData;

    public Reach(IKSegment[] segments, int n)
    {
        minData = new MinReachData(segments, n);
        maxData = new MaxReachData(segments, n);
    }

    public Interval Get(IKSegment[] segments, float angle)
    {
        return new Interval(minData.Get(segments, angle), maxData.Get(segments, angle));
    }
    public Interval ValidInterval { get { return maxData.ValidInterval; } }
}
[System.Serializable]
public class MinReachData
{
    private float[] angles;
    private MinReachRegion[] regions;
    private int n;

    public MinReachData(IKSegment[] segments, int n)
    {
        this.n = n;

        if (n == 1)
        {
            this.angles = new float[2] { segments[0].angleLimit.min, segments[0].angleLimit.max };
            this.regions = new MinReachRegion[1] { new MinReachRegion(0, 1, 0) };
            return;
        }

        float angleSum;
        List<float> angles = new List<float>();
        List<MinReachRegion> regions = new List<MinReachRegion>();
        // [0]: shorter bend side
        // [1]: longer bend side
        int shorterIndex = 0;
        Vector3[] bendPoints = new Vector3[2] 
        { 
            IKUtility.FollowChain(segments, n, IntervalSide.Min),
            IKUtility.FollowChain(segments, n, IntervalSide.Max) 
        };
        float[] lengths = new float[2] 
        { 
            bendPoints[0].magnitude, 
            bendPoints[1].magnitude 
        };
        if (lengths[0] > lengths[1])
        {
            Vector3 temp2 = bendPoints[0];
            bendPoints[0] = bendPoints[1];
            bendPoints[1] = temp2;
            float temp = lengths[0];
            lengths[0] = lengths[1];
            lengths[1] = temp;
            shorterIndex = 1;
        }

        // add angle to full bend on shorter side
        angles.Add(IKUtility.ToAngle(bendPoints[shorterIndex]));
        //regions.Add(new MinReachInterval(shorterIndex, 0));

        Vector3[] prevBends = new Vector3[2] { bendPoints[0], bendPoints[1] };
        Vector3[] lastBends = new Vector3[2] { bendPoints[0], bendPoints[1] };
        int[] k = new int[2] { 1, 1 };
        // start reverse bending by one chain, starting with shorter side
        while (true)
        {
            int activeSide = shorterIndex;
            int inactiveSide = (shorterIndex + 1) % 2;
            int activeK = k[activeSide];

            // calculate the bend position (Forward Kinematics)
            angleSum = 0;
            Vector2 reverseBend = Vector2.zero;
            for (int i = 0; i < activeK; i++)
            {
                angleSum += segments[i].angleLimit[inactiveSide];
                reverseBend += IKUtility.ToVector(angleSum, segments[i].length);
            }
            for (int i = activeK; i < n; i++)
            {
                angleSum += segments[i].angleLimit[activeSide];
                reverseBend += IKUtility.ToVector(angleSum, segments[i].length);
            }

            // update bend history
            prevBends[activeSide] = lastBends[activeSide];
            lastBends[activeSide] = reverseBend;
            // arc angle interval from prevBend to lastBend
            Interval activeInterval = IKUtility.AngleInterval(prevBends[activeSide], lastBends[activeSide]);
            Interval inactiveInterval = IKUtility.AngleInterval(prevBends[inactiveSide], lastBends[inactiveSide]);
            // radius of above interval at end point that is more bent
            float activeLength = lastBends[activeSide].magnitude;
            float inactiveLength = lastBends[inactiveSide].magnitude;

            if (activeLength > inactiveLength)
            {
                break;
            }
            // continue bending the same side
            else
            {
                angles.Add(activeInterval[inactiveSide]);
                regions.Add(new MinReachRegion(inactiveSide, activeSide, activeK - 1));
                k[activeSide]++;
            }

            // infinite loop protection while still writing this function
            if (k[0] > 100 || k[1] > 100)
            {
                break;
            }
        }

        int anglesInsertIndex = angles.Count;

        while (true)
        {
            int activeSide = (shorterIndex + 1) % 2;
            int inactiveSide = shorterIndex;
            int activeK = k[activeSide];

            // calculate the bend position (Forward Kinematics)
            angleSum = 0;
            Vector2 reverseBend = Vector2.zero;
            for (int i = 0; i < activeK; i++)
            {
                angleSum += segments[i].angleLimit[inactiveSide];
                reverseBend += IKUtility.ToVector(angleSum, segments[i].length);
            }
            for (int i = activeK; i < n; i++)
            {
                angleSum += segments[i].angleLimit[activeSide];
                reverseBend += IKUtility.ToVector(angleSum, segments[i].length);
            }

            // update bend history
            prevBends[activeSide] = lastBends[activeSide];
            lastBends[activeSide] = reverseBend;
            // arc angle interval from prevBend to lastBend
            Interval activeInterval = IKUtility.AngleInterval(prevBends[activeSide], lastBends[activeSide]);
            Interval inactiveInterval = IKUtility.AngleInterval(prevBends[inactiveSide], lastBends[inactiveSide]);

            if (Interval.AreOverlapping(activeInterval, inactiveInterval))
            {
                angleSum = 0;
                Vector2 c1 = Vector2.zero;
                for (int i = 0; i < activeK; i++)
                {
                    angleSum += segments[i].angleLimit[inactiveSide];
                    c1 += IKUtility.ToVector(angleSum, segments[i].length);
                }
                angleSum = 0;
                Vector2 c2 = Vector2.zero;
                for (int i = 0; i < k[inactiveSide] - 1; i++)
                {
                    angleSum += segments[i].angleLimit[activeSide];
                    c2 += IKUtility.ToVector(angleSum, segments[i].length);
                }
                angleSum = 0;
                Vector2 r1 = Vector2.zero;
                for (int i = activeK; i < n; i++)
                {
                    angleSum += segments[i].angleLimit[activeSide];
                    r1 += IKUtility.ToVector(angleSum, segments[i].length);
                }
                angleSum = 0;
                Vector2 r2 = Vector2.zero;
                for (int i = k[inactiveSide] - 1; i < n; i++)
                {
                    angleSum += segments[i].angleLimit[inactiveSide];
                    r2 += IKUtility.ToVector(angleSum, segments[i].length);
                }
                //Debug.DrawCircle(c1, r1.magnitude, 128, Color.white, 10);
                //Debug.DrawCircle(c2, r2.magnitude, 128, Color.black, 10);
                Vector2 intersection = IKUtility.TwoCircleIntersection(c1, c2, r1.magnitude, r2.magnitude);

                if (IKUtility.AngleInterval(prevBends[inactiveSide], lastBends[activeSide]).Contains(IKUtility.ToAngle(intersection)))
                {
                    angles.Insert(anglesInsertIndex, IKUtility.ToAngle(lastBends[activeSide]));
                    regions.Insert(anglesInsertIndex - 1, new MinReachRegion(inactiveSide, activeSide, activeK));
                    angles.Insert(anglesInsertIndex, IKUtility.ToAngle(intersection));
                    regions.Insert(anglesInsertIndex - 1, new MinReachRegion(activeSide, inactiveSide, k[inactiveSide] - 1));
                }
                else
                {
                    // TODO:
                }

                // add angle to full bend on longer side
                angles.Add(IKUtility.ToAngle(bendPoints[(shorterIndex + 1) % 2]));
                regions.Add(new MinReachRegion(shorterIndex, (shorterIndex + 1) % 2, 0));
                break;
            }
            else
            {
                angles.Insert(anglesInsertIndex, activeInterval[inactiveSide]);
                regions.Insert(anglesInsertIndex - 1, new MinReachRegion(inactiveSide, activeSide, activeK - 1));
                k[activeSide]++;
            }

            // infinite loop protection while still writing this function
            if (k[0] > 100 || k[1] > 100)
            {
                break;
            }
        }

        // save calculated values
        this.angles = angles.ToArray(); ;
        this.regions = regions.ToArray();
    }

    public float Get(IKSegment[] segments, float angle)
    {
        for (int i = 0; i < regions.Length; i++)
        {
            if (new Interval(angles[i], angles[i + 1]).Contains(angle))
            {
                float angleSum = 0;
                Vector2 center = Vector2.zero;
                for (int j = 0; j < regions[i].MoveChainIndex; j++)
                {
                    angleSum += segments[j].angleLimit[regions[i].StartSide];
                    center += IKUtility.ToVector(angleSum, segments[j].length);
                }
                angleSum = 0;
                Vector2 movingPart = Vector2.zero;
                for (int j = regions[i].MoveChainIndex; j < n; j++)
                {
                    angleSum += segments[j].angleLimit[regions[i].ContinueSide];
                    movingPart += IKUtility.ToVector(angleSum, segments[j].length);
                }

                return IKUtility.IntersectionLineCircle(angle, center, movingPart.magnitude).magnitude;
            }
        }
        return float.NaN;
    }

    public Interval ValidInterval { get { return new Interval(angles[0], angles[angles.Length - 1]); } }
}
[System.Serializable]
public class MaxReachData
{
    private Interval[] angleSums;

    public MaxReachData(IKSegment[] segments, int n)
    {
        angleSums = new Interval[n];
        for (int i = 0; i < n; i++)
        {
            angleSums[i] = IKUtility.AngleInterval(
                IKUtility.FollowChain(segments, i, n, IntervalSide.Min), 
                IKUtility.FollowChain(segments, i, n, IntervalSide.Max));
        }
    }

    public float Get(IKSegment[] segments, float angle)
    {
        float[] lengths = new float[angleSums.Length];
        for (int i = 0; i < lengths.Length; i++)
        {
            lengths[i] = segments[i].length;
        }
        return Get(segments, angle, angleSums.Length, lengths);
    }
    private float Get(IKSegment[] segments, float angle, int n, float[] lengths)
    {
        // to prevent infinite loops (just in case though it shouldn't happen)
        if (n <= 0) { return 0; }
        else if (n == 1) { return lengths[0]; }
        if (angleSums[n - 2].Contains(angle))
        {
            float[] newLengths = new float[lengths.Length - 1];
            for (int i = 0; i < newLengths.Length; i++)
            {
                newLengths[i] = lengths[i];
            }
            newLengths[newLengths.Length - 1] += lengths[newLengths.Length];
            return Get(segments, angle, n - 1, newLengths);
        }
        else if (new Interval(angleSums[n - 1].min, angleSums[n - 2].min).Contains(angle))
        {
            return IKUtility.IntersectionLineCircle(angle, IKUtility.FollowChain(segments, n - 1, IntervalSide.Min), lengths[n - 1]).magnitude;
        }
        else if (new Interval(angleSums[n - 2].max, angleSums[n - 1].max).Contains(angle))
        {
            return IKUtility.IntersectionLineCircle(angle, IKUtility.FollowChain(segments, n - 1, IntervalSide.Max), lengths[n - 1]).magnitude;
        }
        else
        {
            return float.NaN;
        }
    }

    public Interval ValidInterval { get { return angleSums[angleSums.Length - 1]; } }
}
[System.Serializable]
public struct MinReachRegion
{
    private int startSide;
    private int continueSide;
    private int moveChainIndex;

    public MinReachRegion(int startSide, int continueSide, int moveChainIndex)
    {
        this.startSide = startSide;
        this.continueSide = continueSide;
        this.moveChainIndex = moveChainIndex;
    }

    public int StartSide { get { return startSide; } }
    public int ContinueSide { get { return continueSide; } }
    public int MoveChainIndex { get { return moveChainIndex; } }
}
[System.Serializable]
public class IKSegment
{
    public float length;
    public Interval angleLimit;

    public IKSegment(float length, float angleMin, float angleMax)
    {
        this.length = length;
        this.angleLimit = new Interval(angleMin, angleMax);
    }
}
[System.Serializable]
public struct Interval
{
    public float min;
    public float max;

    public Interval(float min, float max)
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
        return i1.min == i2.min && i1.max == i2.max || i1.min < i2.min && i2.min < i1.max || i1.min < i2.max && i2.max < i1.max;
    }
    public static float ClampTo(float value, Interval interval)
    {
        return Mathf.Clamp(value, interval.min, interval.max);
    }

    public bool Contains(float angle)
    {
        return min <= angle && angle <= max;
    }

    public override string ToString()
    {
        return $"[{min}, {max}]";
    }

    public float this[int i] { get { return (i % 2 == 0) ? min : max; } }
}
