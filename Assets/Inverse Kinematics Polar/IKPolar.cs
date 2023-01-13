using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif

public class IKPolar : MonoBehaviour
{
    public bool debugChainBend = true;
    public bool debugStartBendPoints = true;
    public bool debugEndBendPoints = true;
    public bool debugTargetLine = true;
    public bool debugMinReach = true;
    public bool debugMaxReach = true;

    public IKSegment[] segments;
    private Interval[] startAngleSums; // m[k](count)/M[k](count)
    private Vector2[] startMinBendPoints; // a_k(count)
    private Vector2[] startMaxBendPoints; // A_k(count)
    private Interval[] endAngleSums; // n[k](count)/N[k](count)
    private Vector2[] endMinBendPoints; // b_k(count)
    private Vector2[] endMaxBendPoints; // B_k(count)
    private Vector2[] centerPoints; // c_k(count)

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

        if (debugStartBendPoints)
        {
            for (int i = 0; i < segments.Length; i++)
            {
                Debug.DrawLine(Vector2.zero, startMinBendPoints[i], Color.magenta, 0.01f);
                Debug.DrawLine(Vector2.zero, startMaxBendPoints[i], Color.magenta, 0.01f);
            }
        }
        if (debugEndBendPoints)
        {
            for (int i = 0; i < segments.Length; i++)
            {
                Debug.DrawLine(Vector2.zero, endMinBendPoints[i], Color.red, 0.01f);
                Debug.DrawLine(Vector2.zero, endMaxBendPoints[i], Color.red, 0.01f);
            }
        }

        mousePos = Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, -Camera.main.transform.position.z));
        float angle = Mathf.Asin(mousePos.y / mousePos.magnitude) * Mathf.Rad2Deg;
        if (debugTargetLine)
        {
            Debug.DrawRay(Vector3.zero, mousePos, Color.yellow, 0.01f);
        }
        reach = new Vector2(MinReach(angle), MaxReach(angle));
        Debug.DrawRay(mousePos.normalized * reach.x, mousePos.normalized * Mathf.Max(reach.y - reach.x, Mathf.Epsilon), Color.red, 0.01f);
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
        startAngleSums = new Interval[segments.Length];
        startMinBendPoints = new Vector2[segments.Length];
        startMaxBendPoints = new Vector2[segments.Length];
        endAngleSums = new Interval[segments.Length];
        endMinBendPoints = new Vector2[segments.Length];
        endMaxBendPoints = new Vector2[segments.Length];
        for (int i = 0; i < segments.Length; i++)
        {
            startMinBendPoints[i] = CalculateStartMinBendPoint(i, segments.Length);
            startMaxBendPoints[i] = CalculateStartMaxBendPoint(i, segments.Length);
            startAngleSums[i] = new Interval(Vector2.SignedAngle(Vector2.right, startMinBendPoints[i]), Vector2.SignedAngle(Vector2.right, startMaxBendPoints[i]));
            /*startAngleSums[i] = new Interval(Mathf.Atan2(startMinBendPoints[i].x, startMinBendPoints[i].y) * Mathf.Rad2Deg,
                Mathf.Atan2(startMaxBendPoints[i].x, startMaxBendPoints[i].y) * Mathf.Rad2Deg);*/
            /*startAngleSums[i] = new Interval(Mathf.Acos(startMinBendPoints[i].x / startMinBendPoints[i].magnitude) * Mathf.Rad2Deg,
                Mathf.Acos(startMaxBendPoints[i].x / startMaxBendPoints[i].magnitude) * Mathf.Rad2Deg);*/

            endMinBendPoints[i] = CalculateEndMinBendPoint(i, segments.Length);
            endMaxBendPoints[i] = CalculateEndMaxBendPoint(i, segments.Length);
            endAngleSums[i] = new Interval(Vector2.SignedAngle(Vector2.right, endMinBendPoints[i]), Vector2.SignedAngle(Vector2.right, endMaxBendPoints[i]));
            /*endAngleSums[i] = new Interval(Mathf.Atan2(endMinBendPoints[i].x, endMinBendPoints[i].y) * Mathf.Rad2Deg,
                Mathf.Atan2(endMaxBendPoints[i].x, endMaxBendPoints[i].y) * Mathf.Rad2Deg);*/
            /*endAngleSums[i] = new Interval(Mathf.Acos(endMinBendPoints[i].x / endMinBendPoints[i].magnitude) * Mathf.Rad2Deg,
                Mathf.Acos(endMaxBendPoints[i].x / endMaxBendPoints[i].magnitude) * Mathf.Rad2Deg);*/
        }
        centerPoints = new Vector2[segments.Length];
        centerPoints[0] = Vector2.zero;
        for (int i = 0; i < segments.Length - 1; i++)
        {
            centerPoints[i + 1] = new Vector2(centerPoints[i].x + segments[i].length, 0);
        }
    }
    private Vector2 CalculateStartMinBendPoint(int k, int n) { return CalculateStartBendPoint(k, n, 0); }
    private Vector2 CalculateStartMaxBendPoint(int k, int n) { return CalculateStartBendPoint(k, n, 1); }
    private Vector2 CalculateStartBendPoint(int k, int n, int limitIndex)
    {
        Vector2 p = Vector2.zero;
        float angleSum = 0;
        for (int i = 0; i < k; i++)
        {
            angleSum += segments[i].angleLimit[limitIndex];
            p += new Vector2(Mathf.Cos(angleSum * Mathf.Deg2Rad), Mathf.Sin(angleSum * Mathf.Deg2Rad)) * segments[i].length;
        }
        angleSum += segments[k].angleLimit[limitIndex];
        for (int i = k; i < n; i++)
        {
            p += new Vector2(Mathf.Cos(angleSum * Mathf.Deg2Rad), Mathf.Sin(angleSum * Mathf.Deg2Rad)) * segments[i].length;
        }
        return p;
    }
    private Vector2 CalculateEndMinBendPoint(int k, int n) { return CalculateEndBendPoint(k, n, 0); }
    private Vector2 CalculateEndMaxBendPoint(int k, int n) { return CalculateEndBendPoint(k, n, 1); }
    private Vector2 CalculateEndBendPoint(int k, int n, int limitIndex)
    {
        Vector2 p = Vector2.zero;
        for (int i = 0; i < n - k - 1; i++)
        {
            p += new Vector2(segments[i].length, 0);
        }
        float angleSum = 0;
        for (int i = n - k - 1; i < n; i++)
        {
            angleSum += segments[i].angleLimit[limitIndex];
            p += new Vector2(Mathf.Cos(angleSum * Mathf.Deg2Rad), Mathf.Sin(angleSum * Mathf.Deg2Rad)) * segments[i].length;
        }
        return p;
    }
    private Vector2 FollowStartMinChain(int k) { return FollowStartChain(k, 0); }
    private Vector2 FollowStartMaxChain(int k) { return FollowStartChain(k, 1); }
    private Vector2 FollowStartChain(int k, int limitIndex)
    {
        Vector2 p = Vector2.zero;
        float angleSum = 0;
        for (int i = 0; i < k; i++)
        {
            angleSum += segments[i].angleLimit[limitIndex];
            p += new Vector2(Mathf.Cos(angleSum * Mathf.Deg2Rad), Mathf.Sin(angleSum * Mathf.Deg2Rad)) * segments[i].length;
        }
        return p;
    }
    private Vector2 FollowEndMinChain(int k) { return FollowEndChain(k, 0); }
    private Vector2 FollowEndMaxChain(int k) { return FollowEndChain(k, 1); }
    private Vector2 FollowEndChain(int k, int limitIndex)
    {
        Vector2 p = Vector2.zero;
        float angleSum = 0;
        for (int i = segments.Length - k - 1; i < segments.Length; i++)
        {
            angleSum += segments[i].angleLimit[limitIndex];
            p += new Vector2(Mathf.Cos(angleSum * Mathf.Deg2Rad), Mathf.Sin(angleSum * Mathf.Deg2Rad)) * segments[i].length;
        }
        return p;
    }


    public bool Align(Vector3 target)
    {
        return false;
    }
    public float MaxReach(float angle)
    {
        float[] lengths = new float[segments.Length];
        for (int i = 0; i < lengths.Length; i++)
        {
            lengths[i] = segments[i].length;
        }
        return MaxReach(angle, segments.Length, lengths);
    }
    private float MaxReach(float angle, int n, float[] lengths)
    {
        // to prevent infinite loops (just in case though it shouldn't happen)
        if (n <= 0) { return 0; }
        else if (n == 1) { return lengths[0]; }
        if (startAngleSums[n - 2].Contains(angle))
        {
            float[] newLengths = new float[lengths.Length - 1];
            for (int i = 0; i < newLengths.Length; i++)
            {
                newLengths[i] = lengths[i];
            }
            newLengths[newLengths.Length - 1] += lengths[newLengths.Length];
            return MaxReach(angle, n - 1, newLengths);
        }
        else if (new Interval(startAngleSums[n - 1].min, startAngleSums[n - 2].min).Contains(angle))
        {
            Vector3 center = FollowStartMinChain(n - 1);
            if (debugMaxReach)
            {
                Debug.DrawCircle(center, lengths[n - 1], 32, new Color(1.0f, 0.0f, 0.0f));
            }
            return IntersectionLineCircle(angle, center, lengths[n - 1]).magnitude;
        }
        else if (new Interval(startAngleSums[n - 2].max, startAngleSums[n - 1].max).Contains(angle))
        {
            Vector3 center = FollowStartMaxChain(n - 1);
            if (debugMaxReach)
            {
                Debug.DrawCircle(center, lengths[n - 1], 32, new Color(1.0f, 0.25f, 0.0f));
            }
            return IntersectionLineCircle(angle, center, lengths[n - 1]).magnitude;
        }
        else
        {
            return float.NaN;
        }
    }
    public float MinReach(float angle)
    {
        sbyte sign = 0;
        int k = 0;
        Interval prevInterval = new Interval(0, 0);
        for (int n = 0; n < segments.Length; n++)
        {
            if (new Interval(endAngleSums[n].min, prevInterval.min).Contains(angle))
            {
                k = n;
                sign = -1;
            }
            else if (new Interval(prevInterval.max, endAngleSums[n].max).Contains(angle))
            {
                k = n;
                sign = 1;
            }
            prevInterval = endAngleSums[n];
        }
        Debug.Log(sign + " " + k);
        if (sign == 0) { return float.NaN; }

        if (sign < 0)
        {
            Vector2 p = Vector2.zero;
            float angleSum = 0;
            for (int i = segments.Length - k - 1; i < segments.Length; i++)
            {
                angleSum += segments[i].angleLimit[0];
                p += new Vector2(Mathf.Cos(angleSum * Mathf.Deg2Rad), Mathf.Sin(angleSum * Mathf.Deg2Rad)) * segments[i].length;
            }
            if (debugMinReach)
            {
                Debug.DrawCircle(centerPoints[segments.Length - k - 1], p.magnitude, 32, new Color(0.0f, 0.0f, 0.1f));
            }
            return IntersectionLineCircle(angle, centerPoints[segments.Length - k - 1], p.magnitude).magnitude;
        }
        else
        {
            Vector2 p = Vector2.zero;
            float angleSum = 0;
            for (int i = segments.Length - k - 1; i < segments.Length; i++)
            {
                angleSum += segments[i].angleLimit[1];
                p += new Vector2(Mathf.Cos(angleSum * Mathf.Deg2Rad), Mathf.Sin(angleSum * Mathf.Deg2Rad)) * segments[i].length;
            }
            if (debugMinReach)
            {
                Debug.DrawCircle(centerPoints[segments.Length - k - 1], p.magnitude, 32, new Color(0.0f, 0.75f, 0.1f));
            }
            return IntersectionLineCircle(angle, centerPoints[segments.Length - k - 1], p.magnitude).magnitude;
        }
    }

    public Vector2 IntersectionLineCircle(float angle, Vector2 center, float radius)
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
        this.min = min;
        this.max = max;
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
