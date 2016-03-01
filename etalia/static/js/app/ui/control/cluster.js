define(['jquery', 'app/util/utils'], function ($, Utils) {

    var values = {
        0: {label: 'All'},
        1: {label: 'Purple'},
        2: {label: 'Red'},
        3: {label: 'Orange'},
        4: {label: 'White'}
    };

    var Cluster = function (options) {
        this.config = $.extend({
            debug: false,
            element: '#cluster',
            toggle: '#toggle-cluster',
            selection: '#cluster-selection'
        }, options);

        this.$element = $(this.config.element);
        this.$toggle = $(this.config.toggle);
        this.$selection = $(this.config.selection);

        this.selected = null;
    };

    Cluster.prototype.init = function() {
        var that = this,
            $body = $('body');

        this.$toggle.on('click', function(e) {
            Utils.toggleClass(that.$element, 'opened');

            e.preventDefault();
            return false;
        });
        // Close on click out
        $(window).on('click', function(e) {
            if (0 == $(e.target).closest(that.config.element).length) {
                that.$element.removeClass('opened');
            }
        });

        this.$element.on('click', '.choices a', function(e) {
            var value = $(e.target).closest('a').data('cluster');
            that.setValue(value);

            that.$toggle.trigger('click');

            $body.trigger('etalia.control.timespan.change', {
                value: value,
                label: that.getValueColor(value)
            });

            e.preventDefault();
            return false;
        });

        return this;
    };

    Cluster.prototype.getValueColor = function(value) {
        if (values.hasOwnProperty(value)) {
            return values[value].label;
        }
        throw 'Unexpected cluster value';
    };

    Cluster.prototype.setValue = function(value) {
        value = parseInt(value);

        // Deselect current
        this.$selection.removeClass('cluster-' + this.selected);

        // Select 'none/all'
        if (value == 0) {
            this.$selection.hide();
            this.$toggle.find('.cluster-icon').show();
        // Select color
        } else if (0 < value && value <= 4) {
            this.$selection
                .css({display: 'inline-block'})
                .addClass('cluster-' + (this.selected - 1));

            this.$toggle.find('.cluster-icon').hide();
        } else {
            throw 'Unexpected cluster value';
        }

        // Store value
        this.$selection.data('value', value);
        this.selected = value;

        return this;
    };

    Cluster.prototype.getValue = function() {
        return this.$selection.data('value');
    };

    return Cluster;
});
