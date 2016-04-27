define([
    'app',
    'text!app/templates/list.html'
], function (App, template) {

    return App.View.List = App.Backbone.View.extend({
        tagName: 'div',
        class: 'thumb-list',

        template: App.Handlebars.compile(template),

        events: {
            "click #thread-next-page": "onNextPageClick"
        },

        render: function () {

            this.$el.html(this.template({}));

            this.model.view
                .render()
                .appendTo(this.$('.document .inner'));

            return this;
        }
    });
});
