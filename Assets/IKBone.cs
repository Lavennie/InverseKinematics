using UnityEditor;
using UnityEngine;

public class IKBone : MonoBehaviour
{
    public float length;
    public Constraint xLimit = new Constraint(-180.0f, 180.0f);
    public Constraint yLimit = new Constraint(-180.0f, 180.0f);
    public Constraint zLimit = new Constraint(-180.0f, 180.0f);

    public Quaternion GetResting()
    {
        return Quaternion.Euler(xLimit.Center, yLimit.Center, zLimit.Center);
    }

    private void OnDrawGizmos()
    {
        const float radius = 0.3f;
        Vector3 start;
        Vector3 end;

        start = transform.parent.rotation * Quaternion.Euler(xLimit.min, 0, 0) * Vector3.forward;
        end = transform.parent.rotation * Quaternion.Euler(xLimit.max, 0, 0) * Vector3.forward;
        Handles.color = Color.red;
        Handles.DrawWireArc(transform.position, transform.parent.right, start, xLimit.max - xLimit.min, radius);
        Handles.DrawLines(new Vector3[] { 
            transform.position + start * radius, transform.position, 
            transform.position, transform.position + end * radius 
        });

        start = transform.parent.rotation * Quaternion.Euler(0, yLimit.min, 0) * Vector3.forward;
        end = transform.parent.rotation * Quaternion.Euler(0, yLimit.max, 0) * Vector3.forward;
        Handles.color = Color.green;
        Handles.DrawWireArc(transform.position, transform.parent.up, start, yLimit.max - yLimit.min, radius);
        Handles.DrawLines(new Vector3[] {
            transform.position + start * radius, transform.position,
            transform.position, transform.position + end * radius
        });

        start = transform.parent.rotation * Quaternion.Euler(0, 0, zLimit.min) * Vector3.down;
        end = transform.parent.rotation * Quaternion.Euler(0, 0, zLimit.max) * Vector3.down;
        Handles.color = Color.blue;
        Handles.DrawWireArc(transform.position, transform.parent.forward, start, zLimit.max - zLimit.min, radius);
        Handles.DrawLines(new Vector3[] {
            transform.position + start * radius, transform.position,
            transform.position, transform.position + end * radius
        });


        Handles.color = Color.white;
    }
}
[System.Serializable]
public struct Constraint
{
    public float min;
    public float max;

    public Constraint(float min, float max)
    {
        this.min = min;
        this.max = max;
    }

    public float Center { get { return (min + max) / 2.0f; } }
}

public struct AngleRegion
{
    public float min;
    public float max;
    public bool empty;

    public AngleRegion(bool empty)
    {
        this.min = 0.0f;
        this.max = 0.0f;
        this.empty = empty;
    }
    public AngleRegion(float min, float max)
    {
        this.min = min;
        this.max = max;
        this.empty = false;
    }
    public AngleRegion(Constraint constraint)
    {
        this.min = constraint.min;
        this.max = constraint.max;
        this.empty = false;
    }

    public void Intersect(AngleRegion other)
    {
        if (empty || other.empty || max < other.min || other.max < min)
        {
            empty = true;
        }
        else
        {
            min = Mathf.Max(min, other.min);
            max = Mathf.Min(max, other.max);
        }
    }

    public override string ToString()
    {
        return $"min: {min}, max: {max}";
    }

    public AngleRegion Empty { get { return new AngleRegion(true); } }
}
