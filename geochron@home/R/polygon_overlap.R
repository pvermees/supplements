# This set of functions was written by ChatGPT!

# Function to compute the area of a polygon using the Shoelace formula
polygon_area <- function(vertices) {
  n <- nrow(vertices)
  area <- 0
  for (i in 1:(n-1)) {
      area <- area +
          (vertices[i,1] * vertices[i+1,2]) -
          (vertices[i+1,1] * vertices[i,2])
  }
  area <- area +
      (vertices[n,1] * vertices[1,2]) -
      (vertices[1,1] * vertices[n,2])
  return(abs(area) / 2)
}

# Function to check if a point is inside
# a polygon using the ray-casting algorithm
point_in_polygon <- function(point, polygon) {
  n <- nrow(polygon)
  x <- point[1]
  y <- point[2]
  inside <- FALSE
  j <- n
  
  for (i in 1:n) {
    xi <- polygon[i,1]
    yi <- polygon[i,2]
    xj <- polygon[j,1]
    yj <- polygon[j,2]
    
    if (((yi > y) != (yj > y)) &&
        (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
      inside <- !inside
    }
    j <- i
  }
  
  return(inside)
}

# Function to compute intersection points of two line segments
line_intersection <- function(p1, p2, q1, q2) {
  A1 <- p2[2] - p1[2]
  B1 <- p1[1] - p2[1]
  C1 <- A1 * p1[1] + B1 * p1[2]
  
  A2 <- q2[2] - q1[2]
  B2 <- q1[1] - q2[1]
  C2 <- A2 * q1[1] + B2 * q1[2]
  
  det <- A1 * B2 - A2 * B1
  
  if (det == 0) {
    return(NULL)  # Parallel lines
  } else {
    x <- (B2 * C1 - B1 * C2) / det
    y <- (A1 * C2 - A2 * C1) / det
    
    if (min(p1[1], p2[1]) <= x && x <= max(p1[1], p2[1]) &&
        min(p1[2], p2[2]) <= y && y <= max(p1[2], p2[2]) &&
        min(q1[1], q2[1]) <= x && x <= max(q1[1], q2[1]) &&
        min(q1[2], q2[2]) <= y && y <= max(q1[2], q2[2])) {
      return(c(x, y))
    } else {
      return(NULL)
    }
  }
}

# Function to compute the intersection polygon of two polygons
polygon_intersection <- function(poly1, poly2) {
  intersection_points <- matrix(numeric(0), ncol = 2)
  
  # Find intersection points between edges of both polygons
  for (i in 1:nrow(poly1)) {
    p1 <- poly1[i, ]
    p2 <- poly1[ifelse(i == nrow(poly1), 1, i + 1), ]
    
    for (j in 1:nrow(poly2)) {
      q1 <- poly2[j, ]
      q2 <- poly2[ifelse(j == nrow(poly2), 1, j + 1), ]
      
      inter <- line_intersection(p1, p2, q1, q2)
      if (!is.null(inter)) {
        intersection_points <- rbind(intersection_points, inter)
      }
    }
  }
  
  # Add points of poly1 inside poly2
  for (i in 1:nrow(poly1)) {
    if (point_in_polygon(poly1[i, ], poly2)) {
      intersection_points <- rbind(intersection_points, poly1[i, ])
    }
  }
  
  # Add points of poly2 inside poly1
  for (i in 1:nrow(poly2)) {
    if (point_in_polygon(poly2[i, ], poly1)) {
      intersection_points <- rbind(intersection_points, poly2[i, ])
    }
  }
  
  # Remove duplicate points and sort them to form a convex polygon
  if (nrow(intersection_points) > 2) {
    intersection_points <- unique(intersection_points)
    
    # Sort points counterclockwise using centroid
    centroid <- colMeans(intersection_points)
    angles <- atan2(intersection_points[,2] - centroid[2],
                    intersection_points[,1] - centroid[1])
    intersection_points <- intersection_points[order(angles), ]
    
    return(intersection_points)
  } else {
    return(NULL)  # No intersection
  }
}

# Wrapper function to compute the area of intersection
polygon_overlap_area <- function(poly1, poly2) {
  inter_poly <- polygon_intersection(poly1, poly2)
  if (!is.null(inter_poly)) {
    return(polygon_area(inter_poly))
  } else {
    return(0)
  }
}

polygon_test <- function(poly1,poly2){
    poly3 <- polygon_intersection(poly1,poly2)
    xlim <- range(c(poly1[,1],poly2[,1]))
    ylim <- range(c(poly1[,2],poly2[,2]))
    plot(xlim,ylim,type='n')
    polygon(poly1,col='red')
    polygon(poly2,col='blue')
    polygon(poly3,col='green')
}
