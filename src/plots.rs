use anyhow::Result;
use bio::io::fasta;
use plotters::prelude::*;

use crate::gc_content::get_gc_content;

pub fn plot1(record: &fasta::Record, window_size: usize) -> Result<()> {
    let figure = SVGBackend::new("gc-content.svg", (800, 250)).into_drawing_area();

    figure.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&figure)
        .margin(10)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(record.desc().unwrap(), ("sans-serif", 14))
        .build_cartesian_2d(0..(record.seq().len()), 0f64..0.75f64)?;

    // configure labels, axes, etc.
    chart
        .configure_mesh()
        .disable_mesh()
        .x_desc("Mbp")
        .x_labels(40)
        .x_label_formatter(&|x| format!("{}", x / 1000000))
        .y_desc("%GC")
        .y_label_formatter(&|y| format!("{:.0}", y * 100.0))
        .draw()?;

    // plot the gc content with 3 different window sizes, and different colors
    for (window_size_factor, color_idx) in [(1, 0.3), (10, 0.5), (100, 1.0)] {
        let gc_content_iter = get_gc_content(
            record.seq(),
            window_size * window_size_factor,
            (window_size * window_size_factor) / 10,
        );

        chart
            .draw_series(LineSeries::new(gc_content_iter, ygb_color(color_idx)))?
            .label((window_size * window_size_factor / 1000).to_string() + "k")
            .legend(move |(x, y)| {
                Rectangle::new([(x, y - 5), (x + 10, y + 5)], ygb_color(color_idx).filled())
            });
    }

    // configure the legend
    chart
        .configure_series_labels()
        .margin(5)
        .position(SeriesLabelPosition::LowerRight)
        .background_style(WHITE.filled())
        .border_style(&BLACK)
        .draw()?;

    figure.present()?;
    Ok(())
}

pub fn plot2(data: &[(usize, f64)]) -> Result<()> {
    let figure = SVGBackend::new("promotor-gc.svg", (800, 250)).into_drawing_area();

    figure.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&figure)
        .margin(10)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption("Average Promotor Region GC-Content", ("sans-serif", 14))
        .build_cartesian_2d(-950i64..50i64, 0.45f64..0.6f64)?;

    // configure labels, axes, etc.
    chart
        .configure_mesh()
        .disable_mesh()
        .x_desc("bp")
        .y_desc("%GC")
        .y_label_formatter(&|y| format!("{:.0}", y * 100.0))
        .draw()?;

    chart.draw_series(LineSeries::new(
        data.iter().map(|(a, b)| (*a as i64 - 1000, *b)),
        ygb_color(1.0),
    ))?;

    Ok(())
}

pub fn plot3(avg_affinity: impl Iterator<Item = (i32, f64)>) -> Result<()> {
    let figure = SVGBackend::new("promotor-nsome-affinity.svg", (800, 250)).into_drawing_area();

    figure.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&figure)
        .margin(10)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            "Average Promotor Region Nucleosome Affinity",
            ("sans-serif", 14),
        )
        .build_cartesian_2d(-1000i32..100i32, 0.0f64..1.5f64)?;

    // configure labels, axes, etc.
    chart
        .configure_mesh()
        .disable_mesh()
        .x_desc("bp")
        .y_desc("MNase-seq signal")
        .draw()?;

    chart.draw_series(LineSeries::new(avg_affinity, ygb_color(1.0)))?;

    Ok(())
}

pub fn plot4(avg_affinity: Vec<(i32, f64)>, tf_name: &str) -> Result<()> {
    let out_file_name = tf_name.to_string() + ".svg";
    let figure = SVGBackend::new(&out_file_name, (800, 250)).into_drawing_area();

    figure.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&figure)
        .margin(10)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(
            tf_name.to_string() + " Binding Site Nucleosome Affinity",
            ("sans-serif", 14),
        )
        .build_cartesian_2d(-500i32..501i32, 0.0f64..2.0f64)?;

    // configure labels, axes, etc.
    chart
        .configure_mesh()
        .disable_mesh()
        .x_desc("bp")
        .y_desc("MNase-seq signal")
        .draw()?;

    chart.draw_series(LineSeries::new(avg_affinity.into_iter(), ygb_color(1.0)))?;

    Ok(())
}

fn ygb_color(idx: f64) -> plotters::style::RGBColor {
    let color = colorous::YELLOW_GREEN_BLUE.eval_continuous(idx).as_tuple();
    plotters::style::RGBColor {
        0: color.0,
        1: color.1,
        2: color.2,
    }
}
