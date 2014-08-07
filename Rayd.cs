using System;
using UnityEngine;

public struct Rayd
{
	const string defaultNumberFormat = "G5";

	public Vector3d origin;
	internal Vector3d m_direction;

	public Vector3d direction {
		get {
			return m_direction;
		}
		set {
//			if(value.IsZero){
//				throw new System.ArgumentOutOfRangeException("Direction vector cannot be zero");
//			}
			m_direction = value.normalized;
		}
	}

	public Rayd(Ray ray){
		this.origin = (Vector3d)ray.origin;
		this.m_direction = (Vector3d)ray.direction;
	}

	public Rayd(Vector3d origin, Vector3d direction){
		this.origin = origin;
		this.direction = direction;
	}

	public override String ToString(){
		return ToString (defaultNumberFormat);
	}
	public string ToString(string numberFormat){
		return origin.ToString(numberFormat) + "->" + direction.ToString(numberFormat);
	}

	public Vector3d GetPoint (double d)
	{
		return origin + this.m_direction * d;
	}
}


