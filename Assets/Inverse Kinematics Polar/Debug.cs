using UnityEngine;

public class Debug : UnityEngine.Debug
{
    public static void DrawCircle(Circle circle, int segments, Color color)
    {
        DrawCircle(circle.center, circle.radius, segments, color);
    }
    public static void DrawCircle(Circle circle, int segments, Color color, float duration)
    {
        DrawCircle(circle.center, circle.radius, segments, color, duration);
    }
    public static void DrawCircle(Vector3 position, float radius, int segments, Color color)
    {
        DrawCircle(position, radius, segments, color, 0);
    }
    public static void DrawCircle(Vector3 position, float radius, int segments, Color color, float duration)
    {
        // If either radius or number of segments are less or equal to 0, skip drawing
        if (radius <= 0.0f || segments <= 0)
        {
            return;
        }

        // Single segment of the circle covers (360 / number of segments) degrees
        float angleStep = (360.0f / segments);

        // Result is multiplied by Mathf.Deg2Rad constant which transforms degrees to radians
        // which are required by Unity's Mathf class trigonometry methods

        angleStep *= Mathf.Deg2Rad;

        // lineStart and lineEnd variables are declared outside of the following for loop
        Vector3 lineStart = Vector3.zero;
        Vector3 lineEnd = Vector3.zero;

        for (int i = 0; i < segments; i++)
        {
            // Line start is defined as starting angle of the current segment (i)
            lineStart.x = Mathf.Cos(angleStep * i);
            lineStart.y = Mathf.Sin(angleStep * i);

            // Line end is defined by the angle of the next segment (i+1)
            lineEnd.x = Mathf.Cos(angleStep * (i + 1));
            lineEnd.y = Mathf.Sin(angleStep * (i + 1));

            // Results are multiplied so they match the desired radius
            lineStart *= radius;
            lineEnd *= radius;

            // Results are offset by the desired position/origin 
            lineStart += position;
            lineEnd += position;

            // Points are connected using DrawLine method and using the passed color
            DrawLine(lineStart, lineEnd, color, duration);
        }

        DrawPoint(position, color, duration, radius / 10);
    }
    public static void DrawAngle(Vector3 position, float angle, float radius, Color color)
    {
        DrawRay(position, new Vector2(Mathf.Cos(angle * Mathf.Deg2Rad), Mathf.Sin(angle * Mathf.Deg2Rad)) * radius, color);
    }
    public static void DrawAngle(Vector3 position, float angle, float radius, Color color, float duration)
    {
        DrawRay(position, new Vector2(Mathf.Cos(angle * Mathf.Deg2Rad), Mathf.Sin(angle * Mathf.Deg2Rad)) * radius, color, duration);
    }
    public static void DrawAngleInterval(Vector3 position, Interval angleInterval, float radius, Color color)
    {
        DrawAngle(position, angleInterval.min, radius, color);
        DrawAngle(position, angleInterval.max, radius, color);
    }
    public static void DrawAngleInterval(Vector3 position, Interval angleInterval, float radius, Color color, float duration)
    {
        DrawAngle(position, angleInterval.min, radius, color, duration);
        DrawAngle(position, angleInterval.max, radius, color, duration);
    }
    public static void DrawPoint(Vector3 position, Color color, float size = 0.1f)
    {
        DrawLine(position - new Vector3(size, 0, 0), position + new Vector3(size, 0, 0), color);
        DrawLine(position - new Vector3(0, size, 0), position + new Vector3(0, size, 0), color);
    }
    public static void DrawPoint(Vector3 position, Color color, float duration, float size = 0.1f)
    {
        DrawLine(position - new Vector3(size, 0, 0), position + new Vector3(size, 0, 0), color, duration);
        DrawLine(position - new Vector3(0, size, 0), position + new Vector3(0, size, 0), color, duration);
    }
    public static void DrawPoints(CircleIntersection intersection, Color color, float size = 0.1f)
    {
        switch (intersection.Type)
        {
            case CircleIntersection.Variant.Touching:
                DrawPoint(intersection.p1, color, size);
                break;
            case CircleIntersection.Variant.Intersect:
                DrawPoint(intersection.p1, color, size);
                DrawPoint(intersection.p2, color, size);
                break;
        }
    }
    public static void DrawPoints(CircleIntersection intersection, Color color, float duration, float size = 0.1f)
    {
        switch (intersection.Type)
        {
            case CircleIntersection.Variant.Touching:
                DrawPoint(intersection.p1, color, duration, size);
                break;
            case CircleIntersection.Variant.Intersect:
                DrawPoint(intersection.p1, color, duration, size);
                DrawPoint(intersection.p2, color, duration, size);
                break;
        }
    }
}
